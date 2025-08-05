"use client"

import { useEffect, useMemo, useRef, useState } from "react"
import { AgGridReact } from "ag-grid-react"
import { ModuleRegistry, AllCommunityModule } from 'ag-grid-community';
import "ag-grid-community/styles/ag-grid.css"
import "ag-grid-community/styles/ag-theme-alpine.css"
import type { ColDef, GridReadyEvent, GridApi } from "ag-grid-community"
import HtmlCellRenderer from "./html-cell-renderer";
import StructureCellRenderer from "./structure-cell-renderer";
import {
  DropdownMenu,
  DropdownMenuTrigger,
  DropdownMenuContent,
  DropdownMenuCheckboxItem,
} from "@/components/ui/dropdown-menu"
import { Button } from "@/components/ui/button"

ModuleRegistry.registerModules([ AllCommunityModule ]);

interface DataTableProps {
  data: any[]
}

export function DataTable({ data }: DataTableProps) {
  const gridRef = useRef<AgGridReact>(null)
  const [gridApi, setGridApi] = useState<GridApi | null>(null)
  const [columnVisibility, setColumnVisibility] = useState<Record<string, boolean>>({})

  useEffect(() => {
    if (data && data.length > 0) {
      const initialVisibility = Object.keys(data[0]).reduce((acc, key) => {
        acc[key] = true;
        return acc;
      }, {} as Record<string, boolean>);
      setColumnVisibility(initialVisibility);
    }
  }, [data]);

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
      }

      // Handle special characters in header names
      if (key.includes("_gt_") || key.includes("_lt_")) {
        colDef.headerComponent = (params: any) => {
          const headerText = key.replace(/_gt_/g, " &gt; ").replace(/_lt_/g, " &lt; ").replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase());
          return <span>{headerText}</span>;
        };
      } else {
        colDef.headerName = key.replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase());
      }

      // Use a special renderer for HTML content
      if (key.toLowerCase() === 'structure') {
        colDef.cellRenderer = HtmlCellRenderer;
        colDef.autoHeight = true;
        colDef.valueFormatter = params => {
          // Transform base64 string to an img tag
          if (params.value && typeof params.value === 'string' && params.value.length > 100) {
            return `<img src="data:image/svg+xml;base64,${params.value}" alt="Molecule Structure" />`;
          }
          return params.value; // Return original value if not a likely base64 svg
        };
      }
      else if (key.toLowerCase().includes("smile") || key.toLowerCase().includes("nitrostable") || key.toLowerCase().includes("amine")) {
        colDef.cellRenderer = HtmlCellRenderer;
        colDef.autoHeight = true;
      } else {
        // For all other columns, use a formatter to ensure 0 is displayed
        colDef.valueFormatter = params => {
          if (params.value === null || params.value === undefined) {
            return '';
          }
          return String(params.value);
        }
      }
      
      if(key.toLowerCase().includes("smile")) {
        colDef.minWidth = 200
      }
      
      if(key.toLowerCase() === 'structure') {
        colDef.minWidth = 150
      }

      if(key.toLowerCase().includes("nitrostable")) {
        colDef.minWidth = 200
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
  }

  const onGridReady = (params: GridReadyEvent) => {
    setGridApi(params.api)
  }

  useEffect(() => {
    if (gridApi) {
      gridApi.sizeColumnsToFit()
    }
  }, [gridApi, data, columnVisibility]);

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
          Displaying {data.length} rows
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
      <div className="ag-theme-alpine w-full h-[600px] border rounded">
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
          rowHeight={130}
          onGridSizeChanged={(params) => {
            params.api.sizeColumnsToFit()
          }}
          onFirstDataRendered={(params) => {
            params.api.sizeColumnsToFit()
          }}
        />
      </div>
    </div>
  )
}