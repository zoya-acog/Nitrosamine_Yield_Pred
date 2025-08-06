"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Button } from "@/components/ui/button"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Upload, Play, Download, Info } from "lucide-react"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { Badge } from "@/components/ui/badge"
import { DataTable } from "@/components/data-table"
// import { ResultsTable } from "@/components/ResultsTable"

import { FileUpload } from "@/components/file-upload"
import Header from "@/components/header"
import Papa from "papaparse"

export default function ChemicalPipelinePage() {
  const [selectedModel, setSelectedModel] = useState<string>("")
  const [uploadedFile, setUploadedFile] = useState<File | null>(null)
  const [isProcessing, setIsProcessing] = useState(false)
  const [results, setResults] = useState<any[] | null>(null)
  const [error, setError] = useState<string>("")
  const [isClient, setIsClient] = useState(false)

  // Set isClient to true when component mounts
  useEffect(() => {
    setIsClient(true)
  }, [])

  // Debug results whenever they change
  useEffect(() => {
    console.log("=== RESULTS STATE CHANGED ===")
    console.log("Results:", results)
    console.log("Results type:", typeof results)
    console.log("Results is array:", Array.isArray(results))
    console.log("Results length:", results?.length)
    if (results && results.length > 0) {
      console.log("First result item:", results[0])
    }
    console.log("=============================")
  }, [results])

  const handleFileUpload = (file: File) => {
    setUploadedFile(file)
    setError("")
    setResults(null)
    console.log("File uploaded:", file.name)
  }

  const handleProcess = async () => {
    if (!uploadedFile || !selectedModel) {
      setError("Please select a file and model before processing.")
      return
    }

    console.log("=== STARTING PROCESS ===")
    console.log("File:", uploadedFile.name)
    console.log("Model:", selectedModel)

    setIsProcessing(true)
    setError("")
    setResults(null)

    try {
      const formData = new FormData()
      formData.append("file", uploadedFile)
      formData.append("model", selectedModel)

      console.log("Sending request to /api/process")
      
      const response = await fetch("/api/process", {
        method: "POST",
        body: formData,
      })

      console.log("=== API RESPONSE ===")
      console.log("Response status:", response.status)
      console.log("Response ok:", response.ok)
      
      const result = await response.json()
      
      console.log("=== PARSED RESULT ===")
      console.log("Full result object:", result)
      console.log("Result type:", typeof result)
      console.log("Result.success:", result.success)
      console.log("Result.data:", result.data)
      console.log("Result.data type:", typeof result.data)
      console.log("Result.data is array:", Array.isArray(result.data))
      console.log("Result.data length:", result.data?.length)
      
      if (result.data && Array.isArray(result.data) && result.data.length > 0) {
        console.log("First data item:", result.data[0])
        console.log("First data item keys:", Object.keys(result.data[0] || {}))
      }
      console.log("===================")

      if (!response.ok) {
        throw new Error(result.error || `HTTP error! status: ${response.status}`)
      }

      if (result.success && result.data) {
        // A more robust function to parse numeric strings into numbers
        const parseNumerics = (data: any[]) => {
          return data.map(row => {
            const newRow: { [key: string]: any } = {};
            for (const key in row) {
              const value = row[key];
              // Skip parsing for columns that are expected to contain HTML
              if (key.toLowerCase().includes('amine') || key.toLowerCase().includes('structure') || key.toLowerCase().includes('smile')) {
                newRow[key] = value;
              } else if (typeof value === 'string' && !isNaN(parseFloat(value)) && isFinite(value as any)) {
                newRow[key] = parseFloat(value);
              } else {
                newRow[key] = value;
              }
            }
            return newRow;
          });
        };
        
        const parsedData = parseNumerics(result.data);
        setResults(parsedData);

      } else {
        console.error("API response indicates failure:", result)
        throw new Error(result.error || "Processing failed - no data returned")
      }
    } catch (err: any) {
      console.error("=== PROCESS ERROR ===", err)
      setError(err.message || "Error processing file. Please try again.")
    } finally {
      setIsProcessing(false)
      console.log("=== PROCESS COMPLETE ===")
    }
  }

  const downloadResults = () => {
    if (!results) {
      console.log("No results to download")
      return
    }

    console.log("Downloading results:", results)
    
    try {
      const csvContent = Papa.unparse(results)
      const blob = new Blob([csvContent], { type: "text/csv" })
      const url = URL.createObjectURL(blob)
      const a = document.createElement("a")
      a.href = url
      a.download = "processed_results.csv"
      a.click()
      URL.revokeObjectURL(url)
      console.log("Download initiated")
    } catch (error) {
      console.error("Error creating download:", error)
      setError("Error creating download file")
    }
  }

  // Show debug info in UI
  const debugInfo = (
    <div className="mt-4 p-4 bg-gray-100 rounded text-xs">
      <div><strong>Debug Info:</strong></div>
      <div>Results: {results ? `Array with ${results.length} items` : 'null'}</div>
      <div>Results type: {typeof results}</div>
      <div>Is processing: {isProcessing.toString()}</div>
      <div>Error: {error || 'none'}</div>
      <div>File: {uploadedFile?.name || 'none'}</div>
      <div>Model: {selectedModel || 'none'}</div>
    </div>
  )

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      <Header />
      <div className="container mx-auto px-4 py-8 pt-32">
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-8">
          {/* Input Panel */}
          <div className="lg:col-span-1">
            <Card className="h-fit">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Upload className="h-5 w-5" />
                  Input Configuration
                </CardTitle>
                <CardDescription>Upload your CSV file and select the processing model</CardDescription>
              </CardHeader>
              <CardContent className="space-y-6">
                <Tabs defaultValue="file-upload" className="w-full">
                  <TabsList className="grid w-full grid-cols-2">
                    <TabsTrigger value="single-smiles">Single SMILES</TabsTrigger>
                    <TabsTrigger value="file-upload">File Upload</TabsTrigger>
                  </TabsList>

                  <TabsContent value="single-smiles" className="space-y-4">
                    <div>
                      <Label htmlFor="smiles">SMILES</Label>
                      <Input id="smiles" placeholder="Enter SMILES string" className="mt-1" />
                    </div>
                  </TabsContent>

                  <TabsContent value="file-upload" className="space-y-4">
                    <FileUpload onFileUpload={handleFileUpload} />
                    {uploadedFile && (
                      <div className="text-sm text-green-600">✓ File uploaded: {uploadedFile.name}</div>
                    )}
                  </TabsContent>
                </Tabs>

                <div>
                  <Label htmlFor="model">Select Model</Label>
                  <Select value={selectedModel} onValueChange={setSelectedModel}>
                    <SelectTrigger className="mt-1">
                      <SelectValue placeholder="Choose processing model" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="rule-based">
                        <div className="flex items-center gap-2">
                          Rule-based                         
                        </div>
                      </SelectItem>
                      <SelectItem value="gat">
                        <div className="flex items-center gap-2">
                          GAT (Graph Attention)
                        </div>
                      </SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                {error && (
                  <Alert variant="destructive">
                    <AlertDescription>{error}</AlertDescription>
                  </Alert>
                )}

                <Button
                  onClick={handleProcess}
                  disabled={isProcessing || !uploadedFile || !selectedModel}
                  className="w-full"
                >
                  {isProcessing ? (
                    <>
                      <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2" />
                      Processing...
                    </>
                  ) : (
                    <>
                      <Play className="h-4 w-4 mr-2" />
                      Process
                    </>
                  )}
                </Button>

                {results && (
                  <Button onClick={downloadResults} variant="outline" className="w-full bg-transparent">
                    <Download className="h-4 w-4 mr-2" />
                    Download Results
                  </Button>
                )}

                {/* Debug section - remove this in production */}
                {process.env.NODE_ENV === 'development' && debugInfo}
              </CardContent>
            </Card>
          </div>

          {/* Results Panel */}
          <div className="lg:col-span-3">
            <Card className="h-fit">
              <CardHeader>
                <CardTitle>Results</CardTitle>
                {/* <CardDescription>
                  {results
                    ? `Showing ${results.length} processed compounds`
                    : "Results will appear here after processing"}
                </CardDescription> */}
              </CardHeader>
              <CardContent className="min-h-[600px]">
                {results && Array.isArray(results) && results.length > 0 ? (
                  <DataTable data={results} />
                ) : (
                  <div className="flex items-center justify-center h-64 text-gray-500">
                    <div className="text-center">
                      <Info className="h-12 w-12 mx-auto mb-4 opacity-50" />
                      <p>Upload a CSV file and select a model to see results</p>
                      {isProcessing && (
                        <div className="mt-4">
                          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500 mx-auto"></div>
                          <p className="mt-2 text-sm">Processing your file...</p>
                        </div>
                      )}
                    </div>
                  </div>
                )}
              </CardContent>
            </Card>
          </div>
        </div>

        {/* About Section */}
        <Card className="mt-8">
          {/* Keep your existing about section content */}
        </Card>
      </div>
    </div>
  )
}


// 'use client';

// import { useState, useEffect } from 'react';
// import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
// import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
// import { Button } from '@/components/ui/button';
// import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
// import { Input } from '@/components/ui/input';
// import { Label } from '@/components/ui/label';
// import { Upload, Play, Download, Info } from 'lucide-react';
// import { Alert, AlertDescription } from '@/components/ui/alert';
// import { FileUpload } from '@/components/file-upload';
// import Header from '@/components/header';
// // import { ResultsTable } from '../components/ResultsTable';
// import { DataTable } from '@/components/data-table';
// import Papa from 'papaparse';

// export default function ChemicalPipelinePage() {
//   const [selectedModel, setSelectedModel] = useState<string>('');
//   const [uploadedFile, setUploadedFile] = useState<File | null>(null);
//   const [isProcessing, setIsProcessing] = useState(false);
//   const [results, setResults] = useState<any[] | null>(null);
//   const [error, setError] = useState<string>('');
//   const [isClient, setIsClient] = useState(false);

//   // Set isClient to true when component mounts
//   useEffect(() => {
//     setIsClient(true);
//   }, []);

//   // Debug results whenever they change
//   useEffect(() => {
//     console.log('=== RESULTS STATE CHANGED ===');
//     console.log('Results:', results);
//     console.log('Results type:', typeof results);
//     console.log('Results is array:', Array.isArray(results));
//     console.log('Results length:', results?.length);
//     if (results && results.length > 0) {
//       console.log('First result item:', results[0]);
//     }
//     console.log('=============================');
//   }, [results]);

//   const handleFileUpload = (file: File) => {
//     setUploadedFile(file);
//     setError('');
//     setResults(null);
//     console.log('File uploaded:', file.name);
//   };

//   const handleProcess = async () => {
//     if (!uploadedFile || !selectedModel) {
//       setError('Please select a file and model before processing.');
//       return;
//     }

//     console.log('=== STARTING PROCESS ===');
//     console.log('File:', uploadedFile.name);
//     console.log('Model:', selectedModel);

//     setIsProcessing(true);
//     setError('');
//     setResults(null);

//     try {
//       const formData = new FormData();
//       formData.append('file', uploadedFile);
//       formData.append('model', selectedModel);

//       console.log('Sending request to /api/process');

//       const response = await fetch('/api/process', {
//         method: 'POST',
//         body: formData,
//       });

//       console.log('=== API RESPONSE ===');
//       console.log('Response status:', response.status);
//       console.log('Response ok:', response.ok);

//       const result = await response.json();

//       console.log('=== PARSED RESULT ===');
//       console.log('Full result object:', result);
//       console.log('Result type:', typeof result);
//       console.log('Result.success:', result.success);
//       console.log('Result.data:', result.data);
//       console.log('Result.data type:', typeof result.data);
//       console.log('Result.data is array:', Array.isArray(result.data));
//       console.log('Result.data length:', result.data?.length);

//       if (result.data && Array.isArray(result.data) && result.data.length > 0) {
//         console.log('First data item:', result.data[0]);
//         console.log('First data item keys:', Object.keys(result.data[0] || {}));
//       }
//       console.log('===================');

//       if (!response.ok) {
//         throw new Error(result.error || `HTTP error! status: ${response.status}`);
//       }

//       if (result.success && result.data) {
//         // A more robust function to parse numeric strings into numbers
//         const parseNumerics = (data: any[]) => {
//           return data.map(row => {
//             const newRow: { [key: string]: any } = {};
//             for (const key in row) {
//               const value = row[key];
//               // Skip parsing for columns that are expected to contain HTML
//               if (key.toLowerCase().includes('amine') || key.toLowerCase().includes('structure') || key.toLowerCase().includes('smile')) {
//                 newRow[key] = value;
//               } else if (typeof value === 'string' && !isNaN(parseFloat(value)) && isFinite(value as any)) {
//                 newRow[key] = parseFloat(value);
//               } else {
//                 newRow[key] = value;
//               }
//             }
//             return newRow;
//           });
//         };

//         const parsedData = parseNumerics(result.data);
//         setResults(parsedData);
//       } else {
//         console.error('API response indicates failure:', result);
//         throw new Error(result.error || 'Processing failed - no data returned');
//       }
//     } catch (err: any) {
//       console.error('=== PROCESS ERROR ===', err);
//       setError(err.message || 'Error processing file. Please try again.');
//     } finally {
//       setIsProcessing(false);
//       console.log('=== PROCESS COMPLETE ===');
//     }
//   };

//   const downloadResults = () => {
//     if (!results) {
//       console.log('No results to download');
//       return;
//     }

//     console.log('Downloading results:', results);

//     try {
//       const csvContent = Papa.unparse(results);
//       const blob = new Blob([csvContent], { type: 'text/csv' });
//       const url = URL.createObjectURL(blob);
//       const a = document.createElement('a');
//       a.href = url;
//       a.download = 'processed_results.csv';
//       a.click();
//       URL.revokeObjectURL(url);
//       console.log('Download initiated');
//     } catch (error) {
//       console.error('Error creating download:', error);
//       setError('Error creating download file');
//     }
//   };

//   // Convert results to CSV format for ResultsTable
//   const resultsCsv = results && Array.isArray(results) && results.length > 0
//     ? Papa.unparse(results)
//     : '';

//   // Show debug info in UI
//   const debugInfo = (
//     <div className="mt-4 p-4 bg-gray-100 rounded text-xs">
//       <div><strong>Debug Info:</strong></div>
//       <div>Results: {results ? `Array with ${results.length} items` : 'null'}</div>
//       <div>Results type: {typeof results}</div>
//       <div>Is processing: {isProcessing.toString()}</div>
//       <div>Error: {error || 'none'}</div>
//       <div>File: {uploadedFile?.name || 'none'}</div>
//       <div>Model: {selectedModel || 'none'}</div>
//     </div>
//   );

//   return (
//     <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
//       <Header />
//       <div className="container mx-auto px-4 py-8 pt-32">
//         <div className="grid grid-cols-1 lg:grid-cols-4 gap-8">
//           {/* Input Panel */}
//           <div className="lg:col-span-1">
//             <Card className="h-fit">
//               <CardHeader>
//                 <CardTitle className="flex items-center gap-2">
//                   <Upload className="h-5 w-5" />
//                   Input Configuration
//                 </CardTitle>
//                 <CardDescription>Upload your CSV file and select the processing model</CardDescription>
//               </CardHeader>
//               <CardContent className="space-y-6">
//                 <Tabs defaultValue="file-upload" className="w-full">
//                   <TabsList className="grid w-full grid-cols-2">
//                     <TabsTrigger value="single-smiles">Single SMILES</TabsTrigger>
//                     <TabsTrigger value="file-upload">File Upload</TabsTrigger>
//                   </TabsList>

//                   <TabsContent value="single-smiles" className="space-y-4">
//                     <div>
//                       <Label htmlFor="smiles">SMILES</Label>
//                       <Input id="smiles" placeholder="Enter SMILES string" className="mt-1" />
//                     </div>
//                   </TabsContent>

//                   <TabsContent value="file-upload" className="space-y-4">
//                     <FileUpload onFileUpload={handleFileUpload} />
//                     {uploadedFile && (
//                       <div className="text-sm text-green-600">✓ File uploaded: {uploadedFile.name}</div>
//                     )}
//                   </TabsContent>
//                 </Tabs>

//                 <div>
//                   <Label htmlFor="model">Select Model</Label>
//                   <Select value={selectedModel} onValueChange={setSelectedModel}>
//                     <SelectTrigger className="mt-1">
//                       <SelectValue placeholder="Choose processing model" />
//                     </SelectTrigger>
//                     <SelectContent>
//                       <SelectItem value="rule-based">
//                         <div className="flex items-center gap-2">
//                           Rule-based
//                         </div>
//                       </SelectItem>
//                       <SelectItem value="gat">
//                         <div className="flex items-center gap-2">
//                           GAT (Graph Attention)
//                         </div>
//                       </SelectItem>
//                     </SelectContent>
//                   </Select>
//                 </div>

//                 {error && (
//                   <Alert variant="destructive">
//                     <AlertDescription>{error}</AlertDescription>
//                   </Alert>
//                 )}

//                 <Button
//                   onClick={handleProcess}
//                   disabled={isProcessing || !uploadedFile || !selectedModel}
//                   className="w-full"
//                 >
//                   {isProcessing ? (
//                     <>
//                       <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2" />
//                       Processing...
//                     </>
//                   ) : (
//                     <>
//                       <Play className="h-4 w-4 mr-2" />
//                       Process
//                     </>
//                   )}
//                 </Button>

//                 {results && (
//                   <Button onClick={downloadResults} variant="outline" className="w-full bg-transparent">
//                     <Download className="h-4 w-4 mr-2" />
//                     Download Results
//                   </Button>
//                 )}

//                 {/* Debug section - remove this in production */}
//                 {process.env.NODE_ENV === 'development' && debugInfo}
//               </CardContent>
//             </Card>
//           </div>

//           {/* Results Panel */}
//           <div className="lg:col-span-3">
//             <Card className="h-fit">
//               <CardHeader>
//                 <CardTitle>Results</CardTitle>
//                 <CardDescription>
//                   {results
//                     ? `Showing ${results.length} processed compounds`
//                     : 'Results will appear here after processing'}
//                 </CardDescription>
//               </CardHeader>
//               <CardContent className="min-h-[600px]">
//                 {results && Array.isArray(results) && results.length > 0 ? (
//                   <DataTable data={resultsCsv} />
//                 ) : (
//                   <div className="flex items-center justify-center h-64 text-gray-500">
//                     <div className="text-center">
//                       <Info className="h-12 w-12 mx-auto mb-4 opacity-50" />
//                       <p>Upload a CSV file and select a model to see results</p>
//                       {isProcessing && (
//                         <div className="mt-4">
//                           <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500 mx-auto"></div>
//                           <p className="mt-2 text-sm">Processing your file...</p>
//                         </div>
//                       )}
//                     </div>
//                   </div>
//                 )}
//               </CardContent>
//             </Card>
//           </div>
//         </div>

//         {/* About Section */}
//         <Card className="mt-8">
//           {/* Keep your existing about section content */}
//         </Card>
//       </div>
//     </div>
//   );
// }