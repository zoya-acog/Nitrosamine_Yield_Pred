
"use client"

import { ICellRendererParams } from 'ag-grid-community';

const HtmlCellRenderer = (params: ICellRendererParams) => {
  const { value, colDef } = params;
  console.log(`HtmlCellRenderer for ${colDef?.field}: ${value?.slice(0, 50)}...`);

  // Check if value is a base64 string (for Structure or Base64_SVG)
  if (value && typeof value === 'string' && value.length > 100) {
    try {
      // Strip any <img> tag if present (for legacy GAT data)
      let cleanedValue = value;
      if (value.includes('data:image')) {
        const match = value.match(/data:image\/(svg\+xml|png);base64,([^"]+)/);
        if (match) {
          cleanedValue = match[2]; // Extract base64 string
        }
      }

      // Validate base64
      atob(cleanedValue);
      
      // Default to SVG for Base64_SVG and updated GAT
      let imgSrc = `data:image/svg+xml;base64,${cleanedValue}`;
      
      // Fallback to PNG for legacy GAT data
      if (colDef?.field?.toLowerCase() === 'structure' && !cleanedValue.startsWith('PD94bWwg')) {
        imgSrc = `data:image/png;base64,${cleanedValue}`;
      }
      
      return (
        <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
          <img 
            src={imgSrc} 
            alt="Molecule Structure" 
            style={{ maxHeight: '180px', maxWidth: '100%' }} 
            onError={() => console.error(`Failed to load image for ${colDef?.field}`)}
          />
        </div>
      );
    } catch (e) {
      console.error(`Invalid base64 for ${colDef?.field}: ${e}`);
      return <div>Invalid image data</div>;
    }
  }

  // Fallback for non-base64 values (e.g., SMILES, text)
  return <div dangerouslySetInnerHTML={{ __html: value || '' }} />;
};

export default HtmlCellRenderer;