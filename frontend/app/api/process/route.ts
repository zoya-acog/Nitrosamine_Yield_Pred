import { type NextRequest, NextResponse } from "next/server";
import { exec } from "child_process";
import { writeFile, readFile, unlink, mkdir } from "fs/promises";
import { tmpdir } from "os";
import { join } from "path";

export async function POST(request: NextRequest) {
  try {
    const formData = await request.formData();
    const file = formData.get("file") as File;
    const model = formData.get("model") as string;

    if (!file || !model) {
      return NextResponse.json({ error: "File and model are required" }, { status: 400 });
    }

    // 1. Save the uploaded file to a temporary location
    const tempDir = join(tmpdir(), "nitrosamine-yield-pred");
    await mkdir(tempDir, { recursive: true });
    const tempFilePath = join(tempDir, file.name);
    const outputDir = join(tempDir, "output");
    await mkdir(outputDir, { recursive: true });

    const fileBuffer = Buffer.from(await file.arrayBuffer());
    await writeFile(tempFilePath, fileBuffer);

    // 2. Execute the script based on the selected model
    let command: string;
    let outputFilePath: string;

    if (model === "rule-based") {
      outputFilePath = join(outputDir, "new_with_atom_id.csv");
      command = `/home/zoya/Nitrosamine_Yield_Pred/nitro/bin/python /home/zoya/Nitrosamine_Yield_Pred/cli.py rule -i ${tempFilePath} -o ${outputDir} --visualize`;
    } else if (model === "gat") {
      outputFilePath = join(outputDir, "gat_predictions.csv");
      command = `/home/zoya/Nitrosamine_Yield_Pred/nitro/bin/python /home/zoya/Nitrosamine_Yield_Pred/cli.py gat -i ${tempFilePath} -m /home/zoya/Nitrosamine_Yield_Pred/os_duplicate_model.pt -o ${outputDir}`;
    } else {
      // Clean up temp file on error
      await unlink(tempFilePath).catch(console.error);
      return NextResponse.json({ error: "Invalid model selected" }, { status: 400 });
    }

    console.log(`Executing command: ${command}`);
    await new Promise<void>((resolve, reject) => {
      exec(command, (error, stdout, stderr) => {
        if (error) {
          console.error(`exec error: ${error}`);
          console.error(`Script stdout: ${stdout}`);
          console.error(`Script stderr: ${stderr}`);
          // Clean up temp file on error
          unlink(tempFilePath).catch(console.error);
          reject(new Error(`Script execution failed: ${stderr}`));
          return;
        }
        console.log(`Script stdout: ${stdout}`);
        console.error(`Script stderr: ${stderr}`);
        resolve();
      });
    });

    // 3. Read the output
    const resultCsv = await readFile(outputFilePath, "utf-8");

    // 4. Cleanup temporary files
    await unlink(tempFilePath);
    await unlink(outputFilePath);

    // Parse CSV to JSON using proper CSV parsing
    const data = parseCsv(resultCsv);

    console.log(`Parsed ${data.length} rows from CSV`);
    console.log("Sample row:", data[0]);

    return NextResponse.json({ success: true, data });
  } catch (error: any) {
    console.error("Processing error:", error);
    return NextResponse.json({ error: error.message || "Internal server error" }, { status: 500 });
  }
}

// Proper CSV parsing function that handles quoted fields and escaped commas
function parseCsv(csvText: string): Record<string, string>[] {
  const lines = csvText.trim().split('\n');
  if (lines.length === 0) return [];

  const headers = parseCSVLine(lines[0]);
  const data: Record<string, string>[] = [];

  for (let i = 1; i < lines.length; i++) {
    const values = parseCSVLine(lines[i]);
    const row: Record<string, string> = {};
    
    headers.forEach((header, index) => {
      row[header] = values[index] ?? "";
    });
    
    data.push(row);
  }

  return data;
}

// Parse a single CSV line handling quoted fields
function parseCSVLine(line: string): string[] {
  const result: string[] = [];
  let current = '';
  let inQuotes = false;
  let i = 0;

  while (i < line.length) {
    const char = line[i];

    if (char === '"') {
      if (inQuotes && line[i + 1] === '"') {
        // Escaped quote
        current += '"';
        i += 2;
      } else {
        // Toggle quote state
        inQuotes = !inQuotes;
        i++;
      }
    } else if (char === ',' && !inQuotes) {
      // Field separator
      result.push(current.trim());
      current = '';
      i++;
    } else {
      current += char;
      i++;
    }
  }

  // Add the last field
  result.push(current.trim());

  return result;
}