# Use a Node.js base image with Python installed
FROM node:20-slim

# Install Python and system dependencies for RDKit and matplotlib
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    libxrender1 \
    libcairo2 \
    libfreetype6 \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Create and activate a virtual environment
RUN python3 -m venv /app/venv
ENV PATH="/app/venv/bin:$PATH"

# Copy frontend package files
COPY frontend/package.json frontend/pnpm-lock.yaml ./

# Install pnpm and frontend dependencies
RUN npm install -g pnpm && pnpm install

# Copy the rest of the frontend files
COPY frontend/ ./

# Copy backend files
COPY cli.py os_duplicate_model.pt test_inf.csv test_set.csv ./
COPY GAT/ ./GAT/
COPY rule_based/ ./rule_based/

# Install Python dependencies in the virtual environment
COPY rule_based/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir torch pytorch-lightning torch-geometric rdkit

# Ensure permissions for scripts and output directories
RUN mkdir -p /app/output && chmod -R 777 /app/output
RUN chmod -R 755 /app/cli.py /app/GAT /app/rule_based

# Build the Next.js app
RUN pnpm build

# Expose the port for Next.js
EXPOSE 3000

# Start the Next.js app
CMD ["pnpm", "start"]