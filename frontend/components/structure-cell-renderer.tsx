"use client"

import { ICellRendererParams } from 'ag-grid-community';

const StructureCellRenderer = (params: ICellRendererParams) => {
  const { value } = params;

  // Check if the value is a valid-looking Base64 string
  if (value && typeof value === 'string' && value.length > 100) {
    const imgSrc = `data:image/svg+xml;base64,${value}`;
    return (
      <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
        <img 
          src={imgSrc} 
          alt="Molecule Structure" 
          style={{ height: '120px', maxWidth: '100%' }} 
        />
      </div>
    );
  }

  // Return null or a placeholder if there's no image
  return null;
};

export default StructureCellRenderer;
