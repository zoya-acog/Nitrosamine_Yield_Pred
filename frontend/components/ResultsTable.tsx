// 'use client';

// import { useState, useCallback, useEffect, useMemo } from 'react';
// import { AgGridReact } from 'ag-grid-react';
// import { ModuleRegistry, AllCommunityModule } from 'ag-grid-community';
// import type { ColDef, GridApi, GridReadyEvent, ICellRendererParams } from 'ag-grid-community';
// import {
//   DropdownMenu,
//   DropdownMenuTrigger,
//   DropdownMenuContent,
//   DropdownMenuCheckboxItem,
// } from '@/components/ui/dropdown-menu';
// import { Button } from '@/components/ui/button';

// // Register AG Grid modules
// ModuleRegistry.registerModules([AllCommunityModule]);

// interface ResultsTableProps {
//   data: string; // CSV content as a string
// }

// // Image Modal Component
// interface ImageModalProps {
//   isOpen: boolean;
//   imageSrc: string;
//   onClose: () => void;
// }

// function ImageModal({ isOpen, imageSrc, onClose }: ImageModalProps) {
//   if (!isOpen) return null;

//   return (
//     <div
//       className="fixed inset-0 z-50 flex items-center justify-center p-4"
//       style={{ backgroundColor: 'rgba(71, 85, 105, 0.9)' }}
//       onClick={onClose}
//     >
//       <div className="relative max-w-4xl max-h-[90vh] overflow-hidden rounded-lg shadow-2xl">
//         <button
//           onClick={onClose}
//           className="absolute top-2 right-2 z-10 bg-white bg-opacity-80 hover:bg-opacity-100 rounded-full p-2 transition-all duration-200"
//           style={{ color: '#475569' }}
//         >
//           <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
//             <path d="M18 6L6 18M6 6L18 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
//           </svg>
//         </button>
//         <img
//           src={imageSrc}
//           alt="Enlarged view"
//           style={{
//             width: '600px',
//             height: '500px',
//             objectFit: 'contain',
//           }}
//           className="max-w-full max-h-[90vh] object-contain"
//           onClick={(e) => e.stopPropagation()}
//         />
//       </div>
//     </div>
//   );
// }

// // HTML Cell Renderer Component
// function HtmlCellRenderer(params: ICellRendererParams) {
//   const value = params.value;
//   if (!value) return null;

//   // Handle Base64 SVG images
//   if (typeof value === 'string' && value.startsWith('data:image/svg+xml;base64,')) {
//     return <img src={value} alt="Structure" style={{ maxWidth: '100px', maxHeight: '100px' }} />;
//   }

//   // Handle other HTML content or plain text
//   return <div dangerouslySetInnerHTML={{ __html: value }} />;
// }

// export function ResultsTable({ data }: ResultsTableProps) {
//   const [gridApi, setGridApi] = useState<GridApi | null>(null);
//   const [columnVisibility, setColumnVisibility] = useState<Record<string, boolean>>({});
//   const [modalImage, setModalImage] = useState<string>('');
//   const [isModalOpen, setIsModalOpen] = useState(false);

//   // Parse CSV content into rows for AGGrid
//   const parseCsv = (csvContent: string) => {
//     const rows = csvContent.trim().split('\n').map((row) => row.split(','));
//     const headers = rows[0] || ['Name', 'SMILES', 'Docking score (kcal/mol)', 'Ligand efficiency'];
//     const rowData: Record<string, string>[] = rows.slice(1).map((row) =>
//       headers.reduce<Record<string, string>>((obj, header, i) => {
//         obj[header] = row[i] || '';
//         return obj;
//       }, {})
//     );
//     return { headers, rowData };
//   };

//   // Initialize column visibility when data changes
//   useEffect(() => {
//     if (data) {
//       const { headers } = parseCsv(data);
//       const initialVisibility = headers.reduce((acc, key) => {
//         acc[key] = true;
//         return acc;
//       }, {} as Record<string, boolean>);
//       setColumnVisibility(initialVisibility);
//     }
//   }, [data]);

//   // Custom cell renderer for numeric values
//   const numericCellRenderer = (params: any) => {
//     const value = parseFloat(params.value);
//     if (isNaN(value)) return params.value;
//     return value.toFixed(2);
//   };

//   // Custom cell renderer for docking scores with color coding
//   const dockingScoreCellRenderer = (params: any) => {
//     const value = parseFloat(params.value);
//     if (isNaN(value)) return params.value;

//     let colorClass = '';
//     if (value <= -10) colorClass = 'text-green-700 bg-green-50 font-semibold';
//     else if (value <= -8) colorClass = 'text-green-600 bg-green-25 font-medium';
//     else if (value <= -6) colorClass = 'text-yellow-700 bg-yellow-50 font-medium';
//     else if (value <= -4) colorClass = 'text-orange-700 bg-orange-50';
//     else colorClass = 'text-red-700 bg-red-50';

//     return (
//       <span className={`px-2 py-1 rounded-md text-sm ${colorClass}`}>
//         {value.toFixed(2)}
//       </span>
//     );
//   };

//   const handleImageClick = (imageSrc: string) => {
//     setModalImage(imageSrc);
//     setIsModalOpen(true);
//   };

//   const closeModal = () => {
//     setIsModalOpen(false);
//     setModalImage('');
//   };

//   // Enhanced HTML Cell Renderer with image click handler
//   const EnhancedHtmlCellRenderer = (params: any) => {
//     return (
//       <div
//         onClick={(e) => {
//           const target = e.target as HTMLElement;
//           if (target.tagName === 'IMG') {
//             const img = target as HTMLImageElement;
//             handleImageClick(img.src);
//           }
//         }}
//         className="cursor-pointer"
//       >
//         <HtmlCellRenderer {...params} />
//       </div>
//     );
//   };

//   // AGGrid column definitions
//   const columnDefs: ColDef[] = useMemo(() => {
//     if (!data || !parseCsv(data).headers.length) return [];

//     const headers = parseCsv(data).headers;

//     return headers.map((header) => {
//       const baseConfig: ColDef = {
//         headerName: header.replace(/_/g, ' ').replace(/\b\w/g, (l) => l.toUpperCase()),
//         field: header,
//         sortable: true,
//         filter: true,
//         resizable: true,
//         minWidth: 120,
//         hide: !columnVisibility[header],
//       };

//       // Special configuration for specific columns
//       if (header === 'Docking score (kcal/mol)') {
//         return {
//           ...baseConfig,
//           cellRenderer: dockingScoreCellRenderer,
//           sort: 'asc',
//           minWidth: 180,
//           comparator: (valueA: string, valueB: string) => parseFloat(valueA) - parseFloat(valueB),
//         };
//       } else if (header === 'Ligand efficiency') {
//         return {
//           ...baseConfig,
//           cellRenderer: numericCellRenderer,
//           minWidth: 150,
//           flex: 1,
//         };
//       } else if (header.toLowerCase() === 'smiles') {
//         return {
//           ...baseConfig,
//           minWidth: 200,
//           flex: 7,
//           cellClass: 'font-mono text-xs',
//         };
//       } else if (header === 'Name') {
//         return {
//           ...baseConfig,
//           minWidth: 80,
//           cellClass: 'font-medium',
//           pinned: 'left',
//         };
//       } else if (header.toLowerCase() === 'structure' || header.toLowerCase() === 'base64_svg') {
//         return {
//           ...baseConfig,
//           cellRenderer: EnhancedHtmlCellRenderer,
//           autoHeight: true,
//           minWidth: 200,
//         };
//       }
//       return baseConfig;
//     });
//   }, [data, columnVisibility]);

//   // Handle grid ready event
//   const onGridReady = useCallback((params: GridReadyEvent) => {
//     setGridApi(params.api);
//     params.api.sizeColumnsToFit();
//   }, []);

//   // Export to CSV
//   const exportToCsv = () => {
//     if (gridApi) {
//       gridApi.exportDataAsCsv({
//         fileName: 'results_table.csv',
//       });
//     }
//   };

//   // Toggle column visibility
//   const toggleColumnVisibility = (field: string) => {
//     setColumnVisibility((prev) => ({ ...prev, [field]: !prev[field] }));
//   };

//   // Parse the CSV data
//   const { headers, rowData } = parseCsv(data);

//   return (
//     <div className="space-y-4">
//       {/* Column Selection Dropdown */}
//       <div className="flex justify-end">
//         <DropdownMenu>
//           <DropdownMenuTrigger asChild>
//             <Button variant="outline" className="bg-white">
//               Select Columns
//             </Button>
//           </DropdownMenuTrigger>
//           <DropdownMenuContent>
//             {headers.map((header) => (
//               <DropdownMenuCheckboxItem
//                 key={header}
//                 checked={columnVisibility[header]}
//                 onCheckedChange={() => toggleColumnVisibility(header)}
//               >
//                 {header.replace(/_/g, ' ').replace(/\b\w/g, (l) => l.toUpperCase())}
//               </DropdownMenuCheckboxItem>
//             ))}
//           </DropdownMenuContent>
//         </DropdownMenu>
//         <Button onClick={exportToCsv} variant="outline" className="ml-2">
//           Export CSV
//         </Button>
//       </div>

//       {/* AG Grid Table */}
//       <div className="ag-theme-alpine-custom w-full h-[600px] border rounded-xl shadow-lg overflow-hidden">
//         <AgGridReact
//           rowData={rowData}
//           columnDefs={columnDefs}
//           defaultColDef={{
//             flex: 1,
//             minWidth: 100,
//             filter: true,
//             sortable: true,
//             resizable: true,
//             cellClass: 'custom-cell',
//           }}
//           pagination={true}
//           paginationPageSize={50}
//           suppressCellFocus={true}
//           enableCellTextSelection={true}
//           onGridReady={onGridReady}
//           domLayout="normal"
//           animateRows={true}
//           rowSelection="multiple"
//           rowHeight={200}
//         />
//       </div>

//       {/* Image Modal */}
//       <ImageModal isOpen={isModalOpen} imageSrc={modalImage} onClose={closeModal} />

//       {/* Custom Styles */}
//       <style jsx global>{`
//         .ag-theme-alpine-custom {
//           --ag-foreground-color: #1f2937;
//           --ag-background-color: #ffffff;
//           --ag-header-background-color: #f8fafc;
//           --ag-header-foreground-color: #1e293b;
//           --ag-row-hover-background-color: #f0f9ff;
//           --ag-selected-row-background-color: #dbeafe;
//           --ag-border-color: #e2e8f0;
//           --ag-header-height: 48px;
//           --ag-row-height: 200px;
//           --ag-font-size: 14px;
//           --ag-font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
//           --ag-cell-horizontal-padding: 16px;
//           --ag-grid-size: 8px;
//           --ag-row-border-color: #f1f5f9;
//         }

//         .ag-theme-alpine-custom .ag-header-cell {
//           background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
//           border-bottom: 2px solid #e2e8f0;
//           font-weight: 600;
//           color: #1e293b;
//           padding: 0 16px;
//           display: flex;
//           align-items: center;
//           border-right: 1px solid #e2e8f0;
//         }

//         .ag-theme-alpine-custom .ag-header-cell:hover {
//           background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%);
//         }

//         .ag-theme-alpine-custom .ag-header-cell-resize {
//           width: 10px;
//           background: #3b82f6;
//           opacity: 0.4;
//           transition: opacity 0.2s ease;
//         }

//         .ag-theme-alpine-custom .ag-header-cell-resize:hover {
//           opacity: 1;
//           cursor: col-resize;
//         }

//         .ag-theme-alpine-custom .ag-cell {
//           padding: 8px 16px;
//           line-height: 26px;
//           border-bottom: 1px solid #f1f5f9;
//           border-right: 1px solid #e2e8f0;
//           display: flex;
//           align-items: center;
//         }

//         .ag-theme-alpine-custom .ag-row:hover {
//           background-color: #f0f9ff;
//           box-shadow: 0 1px 3px rgba(59, 130, 246, 0.1);
//         }

//         .ag-theme-alpine-custom .ag-row-selected {
//           background-color: #dbeafe !important;
//           border-left: 3px solid #3b82f6;
//         }

//         .ag-theme-alpine-custom .ag-floating-filter {
//           padding: 4px 8px;
//           border: 1px solid #e2e8f0;
//           border-radius: 4px;
//           background-color: #ffffff;
//           box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05);
//           transition: all 0.2s ease;
//         }

//         .ag-theme-alpine-custom .ag-floating-filter:focus-within {
//           border-color: #3b82f6;
//           box-shadow: 0 0 0 2px rgba(59, 130, 246, 0.2);
//         }

//         .ag-theme-alpine-custom .ag-floating-filter-input {
//           font-size: 14px;
//           color: #1f2937;
//           outline: none;
//           padding: 2px 6px;
//           width: 100%;
//         }

//         .ag-theme-alpine-custom .ag-floating-filter-button {
//           background: #f8fafc;
//           border: 1px solid #e2e8f0;
//           border-radius: 4px;
//           padding: 2px 6px;
//           color: #1e293b;
//           transition: all 0.2s ease;
//         }

//         .ag-theme-alpine-custom .ag-floating-filter-button:hover {
//           background: #f0f9ff;
//           border-color: #93c5fd;
//           color: #2563eb;
//         }

//         .ag-theme-alpine-custom .ag-paging-panel {
//           padding: 16px 32px;
//           background: linear-gradient(90deg, #f0f9ff 0%, #e0e7ef 100%);
//           border-top: 2px solid #e2e8f0;
//           border-radius: 0 0 1rem 1rem;
//           box-shadow: 0 -2px 8px 0 rgba(59,130,246,0.04);
//           display: flex;
//           align-items: center;
//           justify-content: center;
//           gap: 16px;
//         }

//         .ag-theme-alpine-custom .ag-paging-button {
//           color: #2563eb;
//           background: #f8fafc;
//           border-radius: 9999px;
//           padding: 0.75rem 1.25rem;
//           margin: 0 6px;
//           font-weight: 700;
//           font-size: 1.15rem;
//           box-shadow: 0 2px 8px rgba(59,130,246,0.08);
//           border: 2px solid #e0e7ef;
//           transition: all 0.18s cubic-bezier(.4,0,.2,1);
//           outline: none;
//           display: flex;
//           align-items: center;
//           gap: 0.5rem;
//         }

//         .ag-theme-alpine-custom .ag-paging-button:hover:not(.ag-disabled),
//         .ag-theme-alpine-custom .ag-paging-button:focus:not(.ag-disabled) {
//           background: #dbeafe;
//           color: #1d4ed8;
//           border-color: #93c5fd;
//           box-shadow: 0 4px 16px 0 rgba(59,130,246,0.13);
//           transform: translateY(-2px) scale(1.06);
//         }

//         .ag-theme-alpine-custom .ag-paging-button.ag-disabled {
//           opacity: 0.5;
//           cursor: not-allowed;
//           background: #f1f5f9;
//           color: #94a3b8;
//           border-color: #e0e7ef;
//           box-shadow: none;
//         }

//         .ag-theme-alpine-custom .ag-paging-number {
//           color: #1e293b;
//           font-weight: 700;
//           font-size: 1.1rem;
//           margin: 0 10px;
//           background: #e0e7ef;
//           border-radius: 0.5rem;
//           padding: 10px 18px;
//           box-shadow: 0 1px 2px rgba(59,130,246,0.04);
//           transition: background 0.2s, color 0.2s;
//         }

//         .ag-theme-alpine-custom .ag-paging-number.ag-paging-number-current {
//           background: #2563eb;
//           color: #fff;
//           box-shadow: 0 2px 8px 0 rgba(59,130,246,0.12);
//           border: 2px solid #2563eb;
//           transform: scale(1.10);
//         }

//         .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar,
//         .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar {
//           width: 8px;
//           height: 8px;
//         }

//         .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-track,
//         .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-track {
//           background: #f1f5f9;
//           border-radius: 4px;
//         }

//         .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-thumb,
//         .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-thumb {
//           background: #cbd5e1;
//           border-radius: 4px;
//           transition: background-color 0.2s ease;
//         }

//         .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-thumb:hover,
//         .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-thumb:hover {
//           background: #94a3b8;
//         }
//       `}</style>
//     </div>
//   );
// }



'use client';

import { useState, useCallback, useEffect, useMemo } from 'react';
import { AgGridReact } from 'ag-grid-react';
import { ModuleRegistry, AllCommunityModule } from 'ag-grid-community';
import type { ColDef, GridApi, GridReadyEvent, ICellRendererParams } from 'ag-grid-community';
import {
  DropdownMenu,
  DropdownMenuTrigger,
  DropdownMenuContent,
  DropdownMenuCheckboxItem,
} from '@/components/ui/dropdown-menu';
import { Button } from '@/components/ui/button';

// Register AG Grid modules
ModuleRegistry.registerModules([AllCommunityModule]);

interface ResultsTableProps {
  data: string; // CSV content as a string
}

// Image Modal Component
interface ImageModalProps {
  isOpen: boolean;
  imageSrc: string;
  onClose: () => void;
}

function ImageModal({ isOpen, imageSrc, onClose }: ImageModalProps) {
  if (!isOpen) return null;

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center p-4"
      style={{ backgroundColor: 'rgba(71, 85, 105, 0.9)' }}
      onClick={onClose}
    >
      <div className="relative max-w-4xl max-h-[90vh] overflow-hidden rounded-lg shadow-2xl">
        <button
          onClick={onClose}
          className="absolute top-2 right-2 z-10 bg-white bg-opacity-80 hover:bg-opacity-100 rounded-full p-2 transition-all duration-200"
          style={{ color: '#475569' }}
        >
          <svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
            <path d="M18 6L6 18M6 6L18 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        </button>
        <img
          src={imageSrc}
          alt="Enlarged view"
          style={{
            width: '600px',
            height: '500px',
            objectFit: 'contain',
          }}
          className="max-w-full max-h-[90vh] object-contain"
          onClick={(e) => e.stopPropagation()}
        />
      </div>
    </div>
  );
}

// HTML Cell Renderer Component
function HtmlCellRenderer(params: ICellRendererParams) {
  const value = params.value;
  if (!value) return null;

  // Handle Base64 SVG images
  if (typeof value === 'string' && value.startsWith('data:image/svg+xml;base64,')) {
    return <img src={value} alt="Structure" style={{ maxWidth: '100px', maxHeight: '100px' }} />;
  }

  // Handle other HTML content or plain text
  return <div dangerouslySetInnerHTML={{ __html: value }} />;
}

export function ResultsTable({ data }: ResultsTableProps) {
  const [gridApi, setGridApi] = useState<GridApi | null>(null);
  const [columnVisibility, setColumnVisibility] = useState<Record<string, boolean>>({});
  const [modalImage, setModalImage] = useState<string>('');
  const [isModalOpen, setIsModalOpen] = useState(false);

  // Parse CSV content into rows for AGGrid
  const parseCsv = (csvContent: string) => {
    const rows = csvContent.trim().split('\n').map((row) => row.split(','));
    const headers = rows[0] || ['Name', 'Recovered SMILES (reactant SMILES)', '2D structure', 'Score', 'Likelihood'];
    const rowData: Record<string, string>[] = rows.slice(1).map((row) =>
      headers.reduce<Record<string, string>>((obj, header, i) => {
        obj[header] = row[i] || '';
        return obj;
      }, {})
    );
    return { headers, rowData };
  };

  // Initialize column visibility when data changes
  useEffect(() => {
    if (data) {
      const { headers } = parseCsv(data);
      const initialVisibility = headers.reduce((acc, key) => {
        acc[key] = ['Name', 'Recovered SMILES (reactant SMILES)', '2D structure', 'Score', 'Likelihood'].includes(key);
        return acc;
      }, {} as Record<string, boolean>);
      setColumnVisibility(initialVisibility);
    }
  }, [data]);

  // Custom cell renderer for numeric values
  const numericCellRenderer = (params: any) => {
    const value = parseFloat(params.value);
    if (isNaN(value)) return params.value;
    return value.toFixed(2);
  };

  // Custom cell renderer for docking scores with color coding (adapted for Score)
  const scoreCellRenderer = (params: any) => {
    const value = parseFloat(params.value);
    if (isNaN(value)) return params.value;

    let colorClass = '';
    if (value <= -10) colorClass = 'text-green-700 bg-green-50 font-semibold';
    else if (value <= -8) colorClass = 'text-green-600 bg-green-25 font-medium';
    else if (value <= -6) colorClass = 'text-yellow-700 bg-yellow-50 font-medium';
    else if (value <= -4) colorClass = 'text-orange-700 bg-orange-50';
    else colorClass = 'text-red-700 bg-red-50';

    return (
      <span className={`px-2 py-1 rounded-md text-sm ${colorClass}`}>
        {value.toFixed(2)}
      </span>
    );
  };

  const handleImageClick = (imageSrc: string) => {
    setModalImage(imageSrc);
    setIsModalOpen(true);
  };

  const closeModal = () => {
    setIsModalOpen(false);
    setModalImage('');
  };

  // Enhanced HTML Cell Renderer with image click handler
  const EnhancedHtmlCellRenderer = (params: any) => {
    return (
      <div
        onClick={(e) => {
          const target = e.target as HTMLElement;
          if (target.tagName === 'IMG') {
            const img = target as HTMLImageElement;
            handleImageClick(img.src);
          }
        }}
        className="cursor-pointer"
      >
        <HtmlCellRenderer {...params} />
      </div>
    );
  };

  // AGGrid column definitions
  const columnDefs: ColDef[] = useMemo(() => {
    if (!data || !parseCsv(data).headers.length) return [];

    const headers = parseCsv(data).headers;

    return headers.map((header) => {
      const baseConfig: ColDef = {
        headerName: header.replace(/_/g, ' ').replace(/\b\w/g, (l) => l.toUpperCase()),
        field: header,
        sortable: true,
        filter: true,
        resizable: true,
        minWidth: 150, // Increased to prevent overlap
        hide: !columnVisibility[header],
        suppressMovable: true, // Prevent manual reordering
      };

      // Special configuration for specific columns
      if (header === 'Score') {
        return {
          ...baseConfig,
          cellRenderer: scoreCellRenderer,
          sort: 'asc',
          minWidth: 180,
          flex: 1,
        };
      } else if (header === 'Likelihood') {
        return {
          ...baseConfig,
          cellRenderer: numericCellRenderer,
          minWidth: 180,
          flex: 1,
        };
      } else if (header === 'Recovered SMILES (reactant SMILES)') {
        return {
          ...baseConfig,
          minWidth: 200,
          flex: 2,
          cellClass: 'font-mono text-xs',
        };
      } else if (header === 'Name') {
        return {
          ...baseConfig,
          minWidth: 150,
          cellClass: 'font-medium',
          pinned: 'left',
        };
      } else if (header === '2D structure') {
        return {
          ...baseConfig,
          cellRenderer: EnhancedHtmlCellRenderer,
          autoHeight: true,
          minWidth: 200,
          flex: 2,
        };
      }
      return baseConfig;
    });
  }, [data, columnVisibility]);

  // Handle grid ready event
  const onGridReady = useCallback((params: GridReadyEvent) => {
    setGridApi(params.api);
    params.api.sizeColumnsToFit(); // Adjusts columns to fit viewport
  }, []);

  // Export to CSV
  const exportToCsv = () => {
    if (gridApi) {
      gridApi.exportDataAsCsv({
        fileName: 'results_table.csv',
      });
    }
  };

  // Toggle column visibility
  const toggleColumnVisibility = (field: string) => {
    setColumnVisibility((prev) => ({ ...prev, [field]: !prev[field] }));
  };

  // Parse the CSV data
  const { headers, rowData } = parseCsv(data);

  return (
    <div className="space-y-4">
      {/* Column Selection Dropdown */}
      <div className="flex justify-end">
        <DropdownMenu>
          <DropdownMenuTrigger asChild>
            <Button variant="outline" className="bg-white">
              Select Columns
            </Button>
          </DropdownMenuTrigger>
          <DropdownMenuContent>
            {headers.map((header) => (
              <DropdownMenuCheckboxItem
                key={header}
                checked={columnVisibility[header]}
                onCheckedChange={() => toggleColumnVisibility(header)}
              >
                {header.replace(/_/g, ' ').replace(/\b\w/g, (l) => l.toUpperCase())}
              </DropdownMenuCheckboxItem>
            ))}
          </DropdownMenuContent>
        </DropdownMenu>
        <Button onClick={exportToCsv} variant="outline" className="ml-2">
          Export CSV
        </Button>
      </div>

      {/* AG Grid Table */}
      <div className="ag-theme-alpine-custom w-full h-[600px] border rounded-xl shadow-lg overflow-hidden">
        <AgGridReact
          rowData={rowData}
          columnDefs={columnDefs}
          defaultColDef={{
            flex: 1,
            minWidth: 150, // Increased to prevent overlap
            filter: true,
            sortable: true,
            resizable: true,
            suppressAutoSizeOnRefresh: true, // Preserve manual resizing
            wrapText: true, // Allow text wrapping
            autoHeight: true, // Adjust row height
            cellClass: 'custom-cell',
          }}
          pagination={true}
          paginationPageSize={50}
          suppressCellFocus={true}
          enableCellTextSelection={true}
          onGridReady={onGridReady}
          domLayout="normal"
          animateRows={true}
          rowSelection="multiple"
          rowHeight={200}
        />
      </div>

      {/* Image Modal */}
      <ImageModal isOpen={isModalOpen} imageSrc={modalImage} onClose={closeModal} />

      {/* Custom Styles */}
      <style jsx global>{`
        .ag-theme-alpine-custom {
          --ag-foreground-color: #1f2937;
          --ag-background-color: #ffffff;
          --ag-header-background-color: #f8fafc;
          --ag-header-foreground-color: #1e293b;
          --ag-row-hover-background-color: #f0f9ff;
          --ag-selected-row-background-color: #dbeafe;
          --ag-border-color: #e2e8f0;
          --ag-header-height: 48px;
          --ag-row-height: 200px;
          --ag-font-size: 14px;
          --ag-font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
          --ag-cell-horizontal-padding: 16px;
          --ag-grid-size: 8px;
          --ag-row-border-color: #f1f5f9;
        }

        .ag-theme-alpine-custom .ag-header-cell {
          background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
          border-bottom: 2px solid #e2e8f0;
          font-weight: 600;
          color: #1e293b;
          padding: 0 16px;
          display: flex;
          align-items: center;
          border-right: 1px solid #e2e8f0;
          min-width: 150px; /* Ensure header width */
        }

        .ag-theme-alpine-custom .ag-header-cell:hover {
          background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%);
        }

        .ag-theme-alpine-custom .ag-header-cell-resize {
          width: 10px;
          background: #3b82f6;
          opacity: 0.4;
          transition: opacity 0.2s ease;
        }

        .ag-theme-alpine-custom .ag-header-cell-resize:hover {
          opacity: 1;
          cursor: col-resize;
        }

        .ag-theme-alpine-custom .ag-cell {
          padding: 8px 16px;
          line-height: 26px;
          border-bottom: 1px solid #f1f5f9;
          border-right: 1px solid #e2e8f0;
          display: flex;
          align-items: center;
          min-width: 150px; /* Ensure cell width */
        }

        .ag-theme-alpine-custom .ag-row:hover {
          background-color: #f0f9ff;
          box-shadow: 0 1px 3px rgba(59, 130, 246, 0.1);
        }

        .ag-theme-alpine-custom .ag-row-selected {
          background-color: #dbeafe !important;
          border-left: 3px solid #3b82f6;
        }

        .ag-theme-alpine-custom .ag-floating-filter {
          padding: 4px 8px;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          background-color: #ffffff;
          box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05);
          transition: all 0.2s ease;
        }

        .ag-theme-alpine-custom .ag-floating-filter:focus-within {
          border-color: #3b82f6;
          box-shadow: 0 0 0 2px rgba(59, 130, 246, 0.2);
        }

        .ag-theme-alpine-custom .ag-floating-filter-input {
          font-size: 14px;
          color: #1f2937;
          outline: none;
          padding: 2px 6px;
          width: 100%;
        }

        .ag-theme-alpine-custom .ag-floating-filter-button {
          background: #f8fafc;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          padding: 2px 6px;
          color: #1e293b;
          transition: all 0.2s ease;
        }

        .ag-theme-alpine-custom .ag-floating-filter-button:hover {
          background: #f0f9ff;
          border-color: #93c5fd;
          color: #2563eb;
        }

        .ag-theme-alpine-custom .ag-paging-panel {
          padding: 16px 32px;
          background: linear-gradient(90deg, #f0f9ff 0%, #e0e7ef 100%);
          border-top: 2px solid #e2e8f0;
          border-radius: 0 0 1rem 1rem;
          box-shadow: 0 -2px 8px 0 rgba(59,130,246,0.04);
          display: flex;
          align-items: center;
          justify-content: center;
          gap: 16px;
        }

        .ag-theme-alpine-custom .ag-paging-button {
          color: #2563eb;
          background: #f8fafc;
          border-radius: 9999px;
          padding: 0.75rem 1.25rem;
          margin: 0 6px;
          font-weight: 700;
          font-size: 1.15rem;
          box-shadow: 0 2px 8px rgba(59,130,246,0.08);
          border: 2px solid #e0e7ef;
          transition: all 0.18s cubic-bezier(.4,0,.2,1);
          outline: none;
          display: flex;
          align-items: center;
          gap: 0.5rem;
        }

        .ag-theme-alpine-custom .ag-paging-button:hover:not(.ag-disabled),
        .ag-theme-alpine-custom .ag-paging-button:focus:not(.ag-disabled) {
          background: #dbeafe;
          color: #1d4ed8;
          border-color: #93c5fd;
          box-shadow: 0 4px 16px 0 rgba(59,130,246,0.13);
          transform: translateY(-2px) scale(1.06);
        }

        .ag-theme-alpine-custom .ag-paging-button.ag-disabled {
          opacity: 0.5;
          cursor: not-allowed;
          background: #f1f5f9;
          color: #94a3b8;
          border-color: #e0e7ef;
          box-shadow: none;
        }

        .ag-theme-alpine-custom .ag-paging-number {
          color: #1e293b;
          font-weight: 700;
          font-size: 1.1rem;
          margin: 0 10px;
          background: #e0e7ef;
          border-radius: 0.5rem;
          padding: 10px 18px;
          box-shadow: 0 1px 2px rgba(59,130,246,0.04);
          transition: background 0.2s, color 0.2s;
        }

        .ag-theme-alpine-custom .ag-paging-number.ag-paging-number-current {
          background: #2563eb;
          color: #fff;
          box-shadow: 0 2px 8px 0 rgba(59,130,246,0.12);
          border: 2px solid #2563eb;
          transform: scale(1.10);
        }

        .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar,
        .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar {
          width: 8px;
          height: 8px;
        }

        .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-track,
        .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-track {
          background: #f1f5f9;
          border-radius: 4px;
        }

        .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-thumb,
        .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-thumb {
          background: #cbd5e1;
          border-radius: 4px;
          transition: background-color 0.2s ease;
        }

        .ag-theme-alpine-custom .ag-body-horizontal-scroll-viewport::-webkit-scrollbar-thumb:hover,
        .ag-theme-alpine-custom .ag-body-vertical-scroll-viewport::-webkit-scrollbar-thumb:hover {
          background: #94a3b8;
        }
      `}</style>
    </div>
  );
}