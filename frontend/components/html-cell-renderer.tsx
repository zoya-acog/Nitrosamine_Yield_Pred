"use client"

import { ICellRendererParams } from 'ag-grid-community';

const HtmlCellRenderer = (params: ICellRendererParams) => {
  return <div dangerouslySetInnerHTML={{ __html: params.value }} />;
};

export default HtmlCellRenderer;
