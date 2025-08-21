import { type NextRequest, NextResponse } from "next/server";
import { exec } from "child_process";
import { writeFile, readFile, unlink, mkdir } from "fs/promises";
import { join } from "path";

export async function POST(request: NextRequest) {
  try {
    const formData = await request.formData();
    const file = formData.get("file") as File;
    const model = formData.get("model") as string;

    if (!file || !model) {
      return NextResponse.json({ error: "File and model are required" }, { status: 400 });
    }

    // 1. Save the uploaded file to a fixed location in /app
    const projectRoot = "/app";
    const outputDir = join(projectRoot, "output");
    await mkdir(outputDir, { recursive: true });
    const tempFilePath = join(outputDir, file.name);

    const fileBuffer = Buffer.from(await file.arrayBuffer());
    await writeFile(tempFilePath, fileBuffer);

    // 2. Execute the script based on the selected model
    let command: string;
    let outputFilePath: string;
    const pythonBin = "/app/venv/bin/python"; // Use virtual environment's Python

    if (model === "rule-based") {
      outputFilePath = join(outputDir, "new_with_atom_id.csv");
      command = `${pythonBin} ${join(projectRoot, "cli.py")} rule -i ${tempFilePath} -o ${outputDir} --visualize`;
    } else if (model === "gat") {
      outputFilePath = join(outputDir, "gat_predictions.csv");
      command = `${pythonBin} ${join(projectRoot, "cli.py")} gat -i ${tempFilePath} -m ${join(projectRoot, "os_duplicate_model.pt")} -o ${outputDir}`;
    } else {
      await unlink(tempFilePath).catch(console.error);
      return NextResponse.json({ error: "Invalid model selected" }, { status: 400 });
    }

    console.log(`Executing command: ${command}`);
    await new Promise<void>((resolve, reject) => {
      exec(command, { cwd: projectRoot }, (error, stdout, stderr) => {
        if (error) {
          console.error(`exec error: ${error}`);
          console.error(`Script stdout: ${stdout}`);
          console.error(`Script stderr: ${stderr}`);
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

    // Parse CSV to JSON
    const data = parseCsv(resultCsv);

    console.log(`Parsed ${data.length} rows from CSV`);
    console.log("Sample row:", data[0]);

    return NextResponse.json({ success: true, data });
  } catch (error: any) {
    console.error("Processing error:", error);
    return NextResponse.json({ error: error.message || "Internal server error" }, { status: 500 });
  }
}

// CSV parsing functions (unchanged)
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

function parseCSVLine(line: string): string[] {
  const result: string[] = [];
  let current = '';
  let inQuotes = false;
  let i = 0;

  while (i < line.length) {
    const char = line[i];
    if (char === '"') {
      if (inQuotes && line[i + 1] === '"') {
        current += '"';
        i += 2;
      } else {
        inQuotes = !inQuotes;
        i++;
      }
    } else if (char === ',' && !inQuotes) {
      result.push(current.trim());
      current = '';
      i++;
    } else {
      current += char;
      i++;
    }
  }
  result.push(current.trim());
  return result;
}