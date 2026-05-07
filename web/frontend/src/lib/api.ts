export type ToolName = 'tm' | 'mc' | 'creator';
export type JobStatus = 'queued' | 'running' | 'succeeded' | 'failed' | 'cancelled';

export type JobRequest = {
  tool: ToolName;
  L: number;
  M: number;
  mode: number;
  span?: number;
  samples?: number;
  run?: number;
  seed?: number;
  options: {
    export_eigenvectors: boolean;
    damping: boolean;
    fugacity?: number;
    section_capacity?: number;
    state_capacity?: number;
    convergence?: number;
  };
};

export type Job = {
  id: string;
  tool: ToolName;
  status: JobStatus;
  command: string[];
  created_at: string;
  started_at?: string;
  finished_at?: string;
  exit_code?: number;
  output_paths: string[];
  log_tail?: string;
  error?: string;
};

export type ResultFile = {
  path: string;
  group: string;
  size_bytes: number;
  modified_at: string;
  meta_path?: string;
  kind: string;
};

export type UofsObject = {
  path: string;
  index: number;
  object_count: number;
  validation_errors: string[];
  polygons: Array<{ start: [number, number, number]; length: number; points: Array<[number, number, number]> }>;
  L: number | null;
  M: number | null;
  span: number | null;
  span_x: number;
  contacts: number;
  shrink_labels: string[];
  linking_number: number | string | null;
};

export type UofsSummary = {
  path: string;
  objects: number;
  polygons: number;
  total_edges: number;
  average_edges_per_object: number;
  average_x_span: number;
  polys_per_object: number;
  L: number | null;
  M: number | null;
  span: number | null;
};

export type DataBrowserItem = {
  name: string;
  path: string;
  is_dir: boolean;
  kind: string;
  size_bytes?: number;
};

export type DataBrowserListing = {
  path: string;
  parent_path?: string;
  items: DataBrowserItem[];
};

export type AnalysisAction = 'summary' | 'validate' | 'count_edges' | 'count_spans' | 'contacts' | 'linking_number' | 'shrink_labels';

export type AnalysisResult = {
  action: AnalysisAction;
  path: string;
  metrics: Record<string, unknown>;
  rows: Array<Record<string, unknown>>;
};

export type TextFile = {
  path: string;
  text: string;
};

async function jsonFetch<T>(url: string, init?: RequestInit): Promise<T> {
  let response: Response;
  try {
    response = await fetch(url, init);
  } catch (err) {
    throw new Error(`Backend unavailable. Start it with "make web-backend" and keep it running. ${err instanceof Error ? err.message : ''}`.trim());
  }
  if (!response.ok) {
    const text = await response.text();
    let detail = text || response.statusText;
    try {
      const parsed = JSON.parse(text);
      detail = parsed.detail ?? detail;
    } catch {
      if (response.status >= 500 && detail.toLowerCase().includes('internal server error')) {
        detail = 'Backend unavailable. Start it with "make web-backend" and keep it running.';
      }
    }
    if (response.status === 404 && url.startsWith('/api/data-browser')) {
      detail = 'Data browser route not found. Restart the app with "make web-dev" so the backend loads the latest code.';
    }
    throw new Error(detail);
  }
  return response.json();
}

export const api = {
  createJob: (request: JobRequest) =>
    jsonFetch<Job>('/api/jobs', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request)
    }),
  jobs: () => jsonFetch<Job[]>('/api/jobs'),
  job: (id: string) => jsonFetch<Job>(`/api/jobs/${id}`),
  cancel: (id: string) => jsonFetch<Job>(`/api/jobs/${id}/cancel`, { method: 'POST' }),
  results: () => jsonFetch<ResultFile[]>('/api/results'),
  browseData: (path: string) => jsonFetch<DataBrowserListing>(`/api/data-browser?path=${encodeURIComponent(path)}`),
  uofsSummary: (path: string) => jsonFetch<UofsSummary>(`/api/uofs/summary?path=${encodeURIComponent(path)}`),
  uofsObject: (path: string, index: number) =>
    jsonFetch<UofsObject>(`/api/uofs/object?path=${encodeURIComponent(path)}&index=${index}`),
  analysis: (path: string, action: AnalysisAction) =>
    jsonFetch<AnalysisResult>(`/api/analysis?path=${encodeURIComponent(path)}&action=${encodeURIComponent(action)}`),
  fileText: (path: string) => jsonFetch<TextFile>(`/api/file-text?path=${encodeURIComponent(path)}`)
};

export function validateRequest(request: JobRequest): string[] {
  const errors: string[] = [];
  if ((request.L + 1) * (request.M + 1) >= 10) {
    errors.push('(L + 1)(M + 1) must be less than 10.');
  }
  if (request.L < 0 || request.M < 0) errors.push('L and M must be non-negative.');
  if (request.mode < 0 || request.mode > 3) errors.push('Mode must be 0, 1, 2, or 3.');
  if ((request.tool === 'mc' || request.tool === 'creator') && (!request.span || request.span < 2 || request.span > 8)) {
    errors.push('Span must be between 2 and 8.');
  }
  if (request.tool === 'mc' && (!request.samples || request.samples < 1 || request.samples > 100)) {
    errors.push('Samples must be between 1 and 100.');
  }
  return errors;
}

export function commandPreview(request: JobRequest): string {
  const bin = request.tool === 'tm' ? 'bin/tm_master' : request.tool === 'mc' ? 'bin/mc_master' : 'bin/creator_all';
  const parts = [bin, '-L', request.L, '-M', request.M, '-m', request.mode];
  if (request.tool !== 'tm') parts.push('-s', request.span ?? 2);
  if (request.tool === 'mc') parts.push('-n', request.samples ?? 1, '-r', request.run ?? 1, '-S', request.seed ?? 1);
  if (request.tool === 'tm') {
    if (request.options.fugacity !== undefined) parts.push('-x', request.options.fugacity);
    if (request.options.export_eigenvectors) parts.push('-E');
    if (request.options.damping) parts.push('-d');
    if (request.options.section_capacity !== undefined) parts.push('-S', request.options.section_capacity);
    if (request.options.state_capacity !== undefined) parts.push('-K', request.options.state_capacity);
    if (request.options.convergence !== undefined) parts.push('-c', request.options.convergence);
  }
  return parts.join(' ');
}
