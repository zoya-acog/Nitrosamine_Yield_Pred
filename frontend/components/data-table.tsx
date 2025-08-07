


"use client"

import { useEffect, useMemo, useRef, useState } from "react"
import { AgGridReact } from "ag-grid-react"
import { ModuleRegistry, AllCommunityModule } from 'ag-grid-community';
import "ag-grid-community/styles/ag-grid.css"
import "ag-grid-community/styles/ag-theme-alpine.css"
import type { ColDef, GridReadyEvent, GridApi } from "ag-grid-community"
import HtmlCellRenderer from "./html-cell-renderer";
import {
  DropdownMenu,
  DropdownMenuTrigger,
  DropdownMenuContent,
  DropdownMenuCheckboxItem,
} from "@/components/ui/dropdown-menu"
import { Button } from "@/components/ui/button"

ModuleRegistry.registerModules([ AllCommunityModule ]);

interface DataTableProps {
  data: any[];
  resultType?: 'rule-based' | 'ml-based' | 'default';
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
            <path d="M18 6L6 18M6 6L18 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
          </svg>
        </button>
        <img 
          src={imageSrc} 
          alt="Enlarged view" 
          style={{ 
            width: '600px',
            height: '500px',
            objectFit: 'contain'
          }}
          className="max-w-full max-h-[90vh] object-contain"
          onClick={(e) => e.stopPropagation()}
        />
      </div>
    </div>
  );
}

export function DataTable({ data, resultType = 'default' }: DataTableProps) {
  const gridRef = useRef<AgGridReact>(null)
  const [gridApi, setGridApi] = useState<GridApi | null>(null)
  const [columnVisibility, setColumnVisibility] = useState<Record<string, boolean>>({})
  const [modalImage, setModalImage] = useState<string>("")
  const [isModalOpen, setIsModalOpen] = useState(false)

  useEffect(() => {
    if (data && data.length > 0) {
      const allKeys = Object.keys(data[0]);
      let initialVisibility: Record<string, boolean>;

      if (resultType === 'rule-based') {
        const defaultVisibleKeys = [
          'name',
          'recovered smiles (reactant smiles)',
          '2d structure',
          'score',
          'likelihood',
          'structure',
          'base64_svg',
        ];
        
        initialVisibility = allKeys.reduce((acc, key) => {
          const lowerKey = key.toLowerCase().replace(/_/g, ' ');
          acc[key] = defaultVisibleKeys.some(defaultKey => defaultKey.includes(lowerKey));
          return acc;
        }, {} as Record<string, boolean>);

      } else {
        initialVisibility = allKeys.reduce((acc, key) => {
          acc[key] = true;
          return acc;
        }, {} as Record<string, boolean>);
      }
      setColumnVisibility(initialVisibility);
    }
  }, [data, resultType]);

  const handleImageClick = (imageSrc: string) => {
    setModalImage(imageSrc);
    setIsModalOpen(true);
  };

  const closeModal = () => {
    setIsModalOpen(false);
    setModalImage("");
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

  const columnDefs: ColDef[] = useMemo(() => {
    if (!data || data.length === 0) {
      return []
    }

    const keys = Object.keys(data[0])

    return keys.map((key) => {
      const colDef: ColDef = {
        field: key,
        sortable: true,
        filter: true,
        resizable: true,
        minWidth: 100,
        hide: !columnVisibility[key],
        tooltipField: key,
        headerTooltip: "Drag to resize column",
      }

      if (key.includes("_gt_") || key.includes("_lt_")) {
        colDef.headerComponent = (params: any) => {
          const headerText = key.replace(/_gt_/g, " &gt; ").replace(/_lt_/g, " &lt; ").replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase());
          return <span>{headerText}</span>;
        };
      } else {
        colDef.headerName = key.replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase());
      }

      if (key.toLowerCase() === 'structure' || key.toLowerCase() === 'base64_svg' || key.toLowerCase() === '2d structure') {
        colDef.cellRenderer = EnhancedHtmlCellRenderer;
        colDef.autoHeight = true;
        colDef.minWidth = 200;
      }
      else if (key.toLowerCase().includes("smile") || key.toLowerCase().includes("nitrostable") || key.toLowerCase().includes("amine")) {
        colDef.cellRenderer = EnhancedHtmlCellRenderer;
        colDef.autoHeight = true;
        colDef.minWidth = 200;
      } else {
        colDef.valueFormatter = params => {
          if (params.value === null || params.value === undefined) {
            return '';
          }
          return String(params.value);
        }
      }

      return colDef
    })
  }, [data, columnVisibility])

  const defaultColDef: ColDef = {
    flex: 1,
    sortable: true,
    filter: true,
    resizable: true,
    floatingFilter: true,
    suppressAutoSize: true, // Preserve manual resizing
  }

  const onGridReady = (params: GridReadyEvent) => {
    setGridApi(params.api)
    params.api.sizeColumnsToFit(); // Initial sizing only
  }

  const exportToCsv = () => {
    if (gridApi) {
      gridApi.exportDataAsCsv({
        fileName: "chemical-analysis-results.csv",
      })
    }
  }

  const toggleColumnVisibility = (field: string) => {
    setColumnVisibility(prev => ({ ...prev, [field]: !prev[field] }));
  };

  if (!data || !Array.isArray(data) || data.length === 0) {
    return (
      <div className="text-center py-8 text-gray-500">
        <div>No data to display</div>
      </div>
    )
  }

  return (
    <div className="space-y-4">
      <div className="flex justify-between items-center mb-2">
        <div className="text-sm text-gray-600">
          Displaying {data.length} molecules
        </div>
        <div className="flex items-center gap-2">
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button variant="outline" className="px-3 py-1 text-sm">
                Select Columns
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent className="max-h-96 overflow-y-auto">
              {Object.keys(columnVisibility).map((field) => (
                <DropdownMenuCheckboxItem
                  key={field}
                  checked={columnVisibility[field]}
                  onCheckedChange={() => toggleColumnVisibility(field)}
                  onSelect={(e) => e.preventDefault()}
                >
                  {field.replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase())}
                </DropdownMenuCheckboxItem>
              ))}
            </DropdownMenuContent>
          </DropdownMenu>
          <Button
            onClick={exportToCsv}
            variant="outline"
            className="px-3 py-1 text-sm"
          >
            Export to CSV
          </Button>
        </div>
      </div>
      <div className="ag-theme-alpine-custom w-full h-[600px] border rounded-xl shadow-lg overflow-hidden">
        <AgGridReact
          ref={gridRef}
          rowData={data}
          columnDefs={columnDefs}
          defaultColDef={defaultColDef}
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
      <ImageModal 
        isOpen={isModalOpen} 
        imageSrc={modalImage} 
        onClose={closeModal} 
      />

      {/* Custom AG Grid Styles */}
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
          border-right: 1px solid #e2e8f0; /* Vertical line for headers */
        }

        .ag-theme-alpine-custom .ag-header-cell:hover {
          background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%);
        }

        /* Resize handle styling */
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
          border-right: 1px solid #e2e8f0; /* Vertical line for cells */
          display: flex;
          align-items: center;
          /* Enable ellipsis for overflowing text */
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
        }

        .ag-theme-alpine-custom .ag-row:hover {
          background-color: #f0f9ff;
          box-shadow: 0 1px 3px rgba(59, 130, 246, 0.1);
        }

        .ag-theme-alpine-custom .ag-row-selected {
          background-color: #dbeafe !important;
          border-left: 3px solid #3b82f6;
        }

        /* Filter styling */
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

        /* Pagination styling */
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

        /* Custom scrollbar */
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
  )
}
