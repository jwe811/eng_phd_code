import React from 'react';
import { createRoot } from 'react-dom/client';
import { Activity, BarChart3, Box, ClipboardCheck, Database, Eye, FileText, Folder, Play, RefreshCw, Search, StopCircle, X } from 'lucide-react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell } from 'recharts';
import { api, AnalysisAction, AnalysisResult, commandPreview, DataBrowserListing, Job, JobRequest, ResultFile, TextFile, UofsSummary, validateRequest } from './lib/api';
import { PolygonPreview } from './components/PolygonPreview';
import './styles.css';

const modes = ['Standard SAP', 'Hamiltonian SAP', '2SAP', '2SAP-Hamiltonian'];
const analysisActions: Array<{ action: AnalysisAction; label: string }> = [
  { action: 'summary', label: 'Summary' },
  { action: 'validate', label: 'Validate' },
  { action: 'count_edges', label: 'Count edges' },
  { action: 'count_spans', label: 'Count spans' },
  { action: 'contacts', label: 'Contacts' },
  { action: 'linking_number', label: 'Linking number' },
  { action: 'shrink_labels', label: 'Shrink labels' }
];

function isSapOr2SapOutput(path: string): boolean {
  const filename = path.split('/').pop() ?? '';
  const lowered = path.toLowerCase();
  if (!path.startsWith('data/') || !filename.endsWith('.txt') || filename.endsWith('.meta')) return false;
  if (lowered.includes('evector') || lowered.includes('transfermatrix')) return false;
  return /^(MCpolys|MC2SAPs|AllSAPs|All2SAPs|AllHamSAPs|AllHam2SAPs)/.test(filename);
}

function isTextFile(path: string): boolean {
  const filename = path.split('/').pop() ?? '';
  return filename.endsWith('.txt') && !filename.endsWith('.meta');
}

function relatedMetaPath(path: string, outputPaths: string[]): string {
  const exact = `${path}.meta`;
  return outputPaths.find((candidate) => candidate === exact) ?? exact;
}

function Skeleton(props: { lines?: number }) {
  return <div className="skeleton-stack">{Array.from({ length: props.lines ?? 3 }, (_, index) => <span key={index} className="skeleton" />)}</div>;
}

function EmptyState(props: { title: string; detail: string }) {
  return <div className="empty-state"><strong>{props.title}</strong><span>{props.detail}</span></div>;
}

function StatusBadge(props: { status: Job['status'] }) {
  return <strong className={`status-badge status-${props.status}`}>{props.status}</strong>;
}

function formatMetric(key: string, value: unknown): string {
  if (key === 'average_edges_per_object' && typeof value === 'number') {
    return value.toFixed(2);
  }
  return String(value);
}

function defaultRequest(tool: JobRequest['tool']): JobRequest {
  return {
    tool,
    L: 1,
    M: 1,
    mode: 0,
    span: tool === 'tm' ? undefined : 2,
    samples: tool === 'mc' ? 1 : undefined,
    run: tool === 'mc' ? 1 : undefined,
    seed: tool === 'mc' ? 1 : undefined,
    options: { export_eigenvectors: false, damping: false }
  };
}

function NumberField(props: { label: string; value?: number; onChange: (value: number | undefined) => void; min?: number; max?: number }) {
  return (
    <label className="field">
      <span>{props.label}</span>
      <input
        type="number"
        value={props.value ?? ''}
        min={props.min}
        max={props.max}
        onChange={(event) => props.onChange(event.currentTarget.value === '' ? undefined : Number(event.currentTarget.value))}
      />
    </label>
  );
}

function Launcher(props: { onJob: (job: Job) => void; onToast: (message: string) => void }) {
  const [tool, setTool] = React.useState<JobRequest['tool']>('creator');
  const [request, setRequest] = React.useState<JobRequest>(defaultRequest('creator'));
  const [error, setError] = React.useState('');
  const [submitting, setSubmitting] = React.useState(false);
  const errors = validateRequest(request);

  function chooseTool(next: JobRequest['tool']) {
    setTool(next);
    setRequest(defaultRequest(next));
    setError('');
  }

  async function submit() {
    setSubmitting(true);
    setError('');
    try {
      const job = await api.createJob(request);
      props.onJob(job);
      props.onToast(`Started ${request.tool} job ${job.id}.`);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setSubmitting(false);
    }
  }

  return (
    <section className="panel launcher">
      <div className="panel-header">
        <div>
          <h2>Run Launcher</h2>
          <p>Guarded local runs with command preview and persisted job history.</p>
        </div>
      </div>
      <div className="tabs">
        {(['creator', 'mc', 'tm'] as const).map((name) => (
          <button key={name} className={tool === name ? 'active' : ''} onClick={() => chooseTool(name)}>
            {name === 'creator' ? 'CreatorAll' : name === 'mc' ? 'Monte Carlo' : 'Transfer Matrix'}
          </button>
        ))}
      </div>
      <div className="grid">
        <NumberField label="L" value={request.L} min={0} onChange={(value) => setRequest({ ...request, L: value ?? 0 })} />
        <NumberField label="M" value={request.M} min={0} onChange={(value) => setRequest({ ...request, M: value ?? 0 })} />
        <label className="field wide">
          <span>Mode</span>
          <select value={request.mode} onChange={(event) => setRequest({ ...request, mode: Number(event.currentTarget.value) })}>
            {modes.map((mode, index) => (
              <option key={mode} value={index}>
                {index}: {mode}
              </option>
            ))}
          </select>
        </label>
        {tool !== 'tm' && <NumberField label="Span" value={request.span} min={2} max={8} onChange={(value) => setRequest({ ...request, span: value })} />}
        {tool === 'mc' && (
          <>
            <NumberField label="Samples" value={request.samples} min={1} max={100} onChange={(value) => setRequest({ ...request, samples: value })} />
            <NumberField label="Run" value={request.run} min={0} onChange={(value) => setRequest({ ...request, run: value })} />
            <NumberField label="Seed" value={request.seed} min={0} onChange={(value) => setRequest({ ...request, seed: value })} />
          </>
        )}
        {tool === 'tm' && (
          <>
            <NumberField label="Fugacity" value={request.options.fugacity} onChange={(value) => setRequest({ ...request, options: { ...request.options, fugacity: value } })} />
            <NumberField label="Section cap" value={request.options.section_capacity} min={1} onChange={(value) => setRequest({ ...request, options: { ...request.options, section_capacity: value } })} />
            <NumberField label="State cap" value={request.options.state_capacity} min={1} onChange={(value) => setRequest({ ...request, options: { ...request.options, state_capacity: value } })} />
            <NumberField label="Convergence" value={request.options.convergence} min={0} onChange={(value) => setRequest({ ...request, options: { ...request.options, convergence: value } })} />
            <label className="check"><input type="checkbox" checked={request.options.export_eigenvectors} onChange={(event) => setRequest({ ...request, options: { ...request.options, export_eigenvectors: event.currentTarget.checked } })} /> Export eigenvectors</label>
            <label className="check"><input type="checkbox" checked={request.options.damping} onChange={(event) => setRequest({ ...request, options: { ...request.options, damping: event.currentTarget.checked } })} /> Damping</label>
          </>
        )}
      </div>
      <div className="command">{commandPreview(request)}</div>
      {errors.length > 0 && <div className="errors">{errors.join(' ')}</div>}
      {error && <div className="errors">{error}</div>}
      <button className="primary" disabled={errors.length > 0 || submitting} onClick={submit}>
        <Play size={16} /> Start run
      </button>
    </section>
  );
}

function JobsPanel(props: { selectedJob?: Job; onSelect: (job: Job) => void; onVisualize: (path: string) => void; onAnalyze: (path: string, action: AnalysisAction) => void; onToast: (message: string) => void }) {
  const [jobs, setJobs] = React.useState<Job[]>([]);
  const [detail, setDetail] = React.useState<Job | undefined>(props.selectedJob);
  const [loading, setLoading] = React.useState(true);
  const [metadata, setMetadata] = React.useState<TextFile | null>(null);

  async function refresh() {
    try {
      const items = await api.jobs();
      setJobs(items);
      const selected = props.selectedJob ?? detail ?? items[0];
      if (selected) setDetail(await api.job(selected.id));
    } finally {
      setLoading(false);
    }
  }

  React.useEffect(() => {
    refresh();
    const timer = window.setInterval(refresh, 2500);
    return () => window.clearInterval(timer);
  }, [props.selectedJob?.id]);

  async function cancel() {
    if (!detail) return;
    const next = await api.cancel(detail.id);
    setDetail(next);
    props.onToast(`Cancelled job ${detail.id}.`);
  }

  async function showMetadata(path: string) {
    try {
      setMetadata(await api.fileText(path));
      props.onToast(`Loaded metadata for ${path}.`);
    } catch (err) {
      props.onToast(err instanceof Error ? err.message : String(err));
    }
  }

  return (
    <section className="panel">
      <div className="panel-header">
        <div>
          <h2>Jobs</h2>
          <p>SQLite-backed local history and logs.</p>
        </div>
      </div>
      <div className="job-list">
        {loading && <Skeleton lines={4} />}
        {!loading && jobs.length === 0 && <EmptyState title="No jobs yet" detail="Start a run from the launcher and it will appear here." />}
        {jobs.map((job) => (
          <button key={job.id} className={detail?.id === job.id ? 'row selected' : 'row'} onClick={async () => { props.onSelect(job); setDetail(await api.job(job.id)); }}>
            <span>{job.tool}</span><StatusBadge status={job.status} /><small>{job.id}</small>
          </button>
        ))}
      </div>
      {detail && (
        <div className="detail">
          <div className="detail-actions">
            <code>{detail.command.join(' ')}</code>
            {detail.status === 'running' && <button className="danger" onClick={cancel}><StopCircle size={16} /> Cancel</button>}
          </div>
          {detail.status === 'succeeded' && detail.output_paths.length > 0 && (
            <div className="quick-actions">
              {detail.output_paths.filter(isSapOr2SapOutput).map((path) => (
                <div className="quick-action-row" key={path}>
                  <span>{path}</span>
                  <button onClick={() => props.onVisualize(path)}>Visualize</button>
                  <button onClick={() => props.onAnalyze(path, 'summary')}>Analysis</button>
                </div>
              ))}
            </div>
          )}
          <pre>{detail.log_tail || 'No logs yet.'}</pre>
        </div>
      )}
    </section>
  );
}

function RawFileViewer(props: { path: string; onClose: () => void }) {
  const [text, setText] = React.useState<string | null>(null);
  const [error, setError] = React.useState('');

  React.useEffect(() => {
    let cancelled = false;
    setText(null);
    setError('');
    api.fileText(props.path)
      .then((result) => { if (!cancelled) setText(result.text); })
      .catch((err) => { if (!cancelled) setError(err instanceof Error ? err.message : String(err)); });
    return () => { cancelled = true; };
  }, [props.path]);

  React.useEffect(() => {
    function onKeyDown(event: KeyboardEvent) {
      if (event.key === 'Escape') props.onClose();
    }
    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  }, [props.onClose]);

  return (
    <div className="raw-viewer-overlay" onClick={(e) => { if (e.target === e.currentTarget) props.onClose(); }}>
      <div className="raw-viewer">
        <div className="raw-viewer-header">
          <div>
            <strong>{props.path.split('/').pop()}</strong>
            <small>{props.path}</small>
          </div>
          <button className="icon" onClick={props.onClose} aria-label="Close viewer"><X size={20} /></button>
        </div>
        {error && <div className="errors">{error}</div>}
        {text === null && !error && <Skeleton lines={12} />}
        {text !== null && (
          <div className="raw-viewer-content">
            <textarea readOnly value={text} spellCheck={false} />
          </div>
        )}
      </div>
    </div>
  );
}

function DataPanel(props: { onVisualize: (path: string) => void; onAnalyze: (path: string, action: AnalysisAction) => void; onToast: (message: string) => void }) {
  const [directory, setDirectory] = React.useState('data');
  const [listing, setListing] = React.useState<DataBrowserListing | null>(null);
  const [loading, setLoading] = React.useState(true);
  const [error, setError] = React.useState('');
  const [viewingFile, setViewingFile] = React.useState<string | null>(null);

  React.useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError('');
    api.browseData(directory)
      .then((next) => {
        if (!cancelled) setListing(next);
      })
      .catch((err) => {
        if (!cancelled) setError(err instanceof Error ? err.message : String(err));
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [directory]);

  const filteredItems = listing?.items.filter(item => item.is_dir || isTextFile(item.path)) || [];

  return (
    <section className="panel">
      <div className="panel-header">
        <div>
          <h2>Data</h2>
          <p>Repo-local SAP and 2SAP files.</p>
        </div>
      </div>
      
      <div className="browser-path-header">
        <button className="icon" onClick={() => setDirectory('data')} title="Go to data root"><Database size={16} /></button>
        <strong>{listing?.path || directory}</strong>
      </div>

      {loading && <Skeleton lines={5} />}
      {error && <div className="errors">{error}</div>}
      
      {!loading && !error && (
        <div className="data-browser-list">
          {listing?.parent_path && (
            <button className="browser-row" onClick={() => setDirectory(listing.parent_path!)}>
              <Folder size={17} /><span>..</span><small>Parent directory</small>
            </button>
          )}
          {filteredItems.length === 0 && <EmptyState title="No files found" detail="No text files were found in this directory." />}
          {filteredItems.map((item) => (
            <div key={item.path} className="data-row">
              <div className="data-info" onClick={() => item.is_dir && setDirectory(item.path)} style={{ cursor: item.is_dir ? 'pointer' : 'default' }}>
                {item.is_dir ? <Folder size={17} /> : <FileText size={17} />}
                <div>
                  <strong>{item.name}</strong>
                  <small>{item.is_dir ? 'Directory' : `${item.kind}${item.size_bytes ? ` · ${item.size_bytes} bytes` : ''}`}</small>
                </div>
              </div>
              {!item.is_dir && (
                <div className="button-strip">
                  <button onClick={() => setViewingFile(item.path)} title="View raw file"><Eye size={15} /></button>
                  {item.kind === 'uofs' && (
                    <>
                      <button onClick={() => props.onVisualize(item.path)}>Visualize</button>
                      <button onClick={() => props.onAnalyze(item.path, 'summary')}>Analysis</button>
                    </>
                  )}
                </div>
              )}
            </div>
          ))}
        </div>
      )}
      {viewingFile && <RawFileViewer path={viewingFile} onClose={() => setViewingFile(null)} />}
    </section>
  );
}

function UofsBrowser(props: { currentPath: string; onChoose: (path: string) => void; onClose: () => void }) {
  const [directory, setDirectory] = React.useState('data');
  const [listing, setListing] = React.useState<DataBrowserListing | null>(null);
  const [error, setError] = React.useState('');

  React.useEffect(() => {
    let cancelled = false;
    setError('');
    api.browseData(directory)
      .then((next) => {
        if (!cancelled) setListing(next);
      })
      .catch((err) => {
        if (!cancelled) setError(err instanceof Error ? err.message : String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [directory]);

  return (
    <div className="browser">
      <div className="browser-header">
        <div>
          <strong>{listing?.path ?? directory}</strong>
          <small>Browse files under data/</small>
        </div>
        <div className="browser-actions">
          <button className="icon" onClick={() => setDirectory('data')} aria-label="Go to data root"><Database size={18} /></button>
          <button className="icon" onClick={props.onClose} aria-label="Close browser"><X size={18} /></button>
        </div>
      </div>
      {error && <div className="errors">{error}</div>}
      <div className="browser-list">
        {listing?.parent_path && (
          <button className="browser-row" onClick={() => setDirectory(listing.parent_path!)}>
            <Folder size={17} /><span>..</span><small>Parent directory</small>
          </button>
        )}
        {listing?.items.map((item) => (
          <button
            key={item.path}
            className={item.path === props.currentPath ? 'browser-row selected' : 'browser-row'}
            onClick={() => {
              if (item.is_dir) {
                setDirectory(item.path);
              } else if (item.kind === 'uofs') {
                props.onChoose(item.path);
              }
            }}
            disabled={!item.is_dir && item.kind !== 'uofs'}
          >
            {item.is_dir ? <Folder size={17} /> : <FileText size={17} />}
            <span>{item.name}</span>
            <small>{item.is_dir ? 'Directory' : `${item.kind}${item.size_bytes ? ` · ${item.size_bytes} bytes` : ''}`}</small>
          </button>
        ))}
      </div>
    </div>
  );
}

function AnalysisChart({ result }: { result: AnalysisResult }) {
  let chartData: any[] = [];
  let xKey = 'value';
  let yKey = 'count';
  let label = 'Count';

  if (result.action === 'count_edges' || result.action === 'count_spans' || result.action === 'contacts') {
    chartData = result.rows;
    label = result.action === 'count_edges' ? 'Edge count' : result.action === 'count_spans' ? 'Span' : 'Contacts';
  } else if (result.action === 'linking_number') {
    const counts: Record<string, number> = {};
    result.rows.forEach((row: any) => {
      const val = row.linking_number;
      if (val !== undefined) {
        counts[val] = (counts[val] || 0) + 1;
      }
    });
    chartData = Object.entries(counts).map(([value, count]) => ({ value, count })).sort((a, b) => Number(a.value) - Number(b.value));
    xKey = 'value';
    label = 'Linking number';
  } else if (result.action === 'shrink_labels') {
    const counts: Record<string, number> = {};
    result.rows.forEach((row: any) => {
      (row.labels || []).forEach((labelStr: string) => {
        counts[labelStr] = (counts[labelStr] || 0) + 1;
      });
    });
    chartData = Object.entries(counts).map(([labelStr, count]) => ({ label: labelStr, count })).sort((a, b) => b.count - a.count);
    xKey = 'label';
    label = 'Label frequency';
  } else {
    return null;
  }

  if (chartData.length === 0) return null;

  return (
    <div className="analysis-chart-container">
      <div className="chart-header">
        <strong>Distribution: {label}</strong>
      </div>
      <div style={{ width: '100%', height: 320 }}>
        <ResponsiveContainer width="100%" height="100%">
          <BarChart data={chartData} margin={{ top: 20, right: 30, left: 0, bottom: 20 }}>
            <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#d9ddcf" />
            <XAxis 
              dataKey={xKey} 
              stroke="#5f665d" 
              fontSize={12} 
              tickLine={false} 
              axisLine={{ stroke: '#d9ddcf' }}
              label={{ value: xKey === 'value' ? label : 'Labels', position: 'insideBottom', offset: -10, fontSize: 11, fill: '#5f665d' }}
            />
            <YAxis 
              stroke="#5f665d" 
              fontSize={12} 
              tickLine={false} 
              axisLine={{ stroke: '#d9ddcf' }}
              label={{ value: 'Frequency', angle: -90, position: 'insideLeft', fontSize: 11, fill: '#5f665d' }}
            />
            <Tooltip 
              cursor={{ fill: 'rgba(11, 106, 65, 0.05)' }}
              contentStyle={{ backgroundColor: '#fff', border: '1px solid #d9ddcf', borderRadius: '4px', boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
              itemStyle={{ color: '#0B6A41', fontWeight: 600 }}
              labelStyle={{ color: '#000', marginBottom: 4, fontWeight: 700 }}
            />
            <Bar dataKey={yKey} fill="#0B6A41" radius={[2, 2, 0, 0]} barSize={40}>
              {chartData.map((_entry, index) => (
                <Cell key={`cell-${index}`} fill={index % 2 === 0 ? '#0B6A41' : '#148c58'} />
              ))}
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}

function AnalysisPanel(props: { path: string; action?: AnalysisAction; onPathChange: (path: string) => void; onToast: (message: string) => void }) {
  const [path, setPath] = React.useState(props.path);
  const [browserOpen, setBrowserOpen] = React.useState(false);
  const [result, setResult] = React.useState<AnalysisResult | null>(null);
  const [activeAction, setActiveAction] = React.useState<AnalysisAction | undefined>(props.action);
  const [loading, setLoading] = React.useState(false);
  const [error, setError] = React.useState('');

  React.useEffect(() => {
    setPath(props.path);
  }, [props.path]);

  React.useEffect(() => {
    if (props.path && props.action) {
      setActiveAction(props.action);
      run(props.action, props.path);
    }
  }, [props.path, props.action]);

  async function run(action: AnalysisAction, targetPath = path) {
    setActiveAction(action);
    setLoading(true);
    setError('');
    setResult(null);
    try {
      const next = await api.analysis(targetPath, action);
      setResult(next);
      props.onToast(`${analysisActions.find((item) => item.action === action)?.label ?? action} complete.`);
    } catch (err) {
      const message = err instanceof Error ? err.message : String(err);
      setError(message);
      props.onToast(message);
    } finally {
      setLoading(false);
    }
  }

  const rows = result?.rows ?? [];
  const columns = rows.length > 0 ? Array.from(new Set(rows.flatMap((row) => Object.keys(row)))) : [];

  return (
    <section className="panel analysis-panel">
      <div className="panel-header">
        <div>
          <h2>Analysis</h2>
          <p>Run postprocessing helpers on UofS files.</p>
        </div>
      </div>
      <label className="field">
        <span>UofS path</span>
        <div className="path-control">
          <input value={path} onChange={(event) => { setPath(event.currentTarget.value); props.onPathChange(event.currentTarget.value); }} />
          <button onClick={() => setBrowserOpen((open) => !open)}><Folder size={16} /> Browse data/</button>
        </div>
      </label>
      {browserOpen && <UofsBrowser currentPath={path} onClose={() => setBrowserOpen(false)} onChoose={(nextPath) => { setPath(nextPath); props.onPathChange(nextPath); setBrowserOpen(false); }} />}
      <div className="analysis-actions">
        {analysisActions.map((item) => (
          <button key={item.action} className={activeAction === item.action ? 'active' : ''} onClick={() => run(item.action)} disabled={loading}>{item.label}</button>
        ))}
      </div>
      {error && <div className="errors">{error}</div>}
      {loading && <Skeleton lines={5} />}
      {!loading && !result && !error && (
        !path ? (
          <EmptyState title="No file selected" detail="Please choose a SAP/2SAP file to analyze." />
        ) : (
          <EmptyState title="No analysis yet" detail="Choose an analysis action above to begin." />
        )
      )}
      {result && !loading && (
        <div className="analysis-output">
          <h3>{analysisActions.find((item) => item.action === result.action)?.label ?? result.action}</h3>
          <div className="metadata-grid">
            {Object.entries(result.metrics)
              .filter(([key]) => key.toLowerCase() !== 'path')
              .map(([key, value]) => (
                <div key={key} className={key.toLowerCase().includes('path') ? 'wide-metric' : undefined}><span>{key}</span><strong>{formatMetric(key, value)}</strong></div>
              ))}
          </div>
          <AnalysisChart result={result} />
          {rows.length > 0 && (
            <div className="table-wrap">
              <table>
                <thead><tr>{columns.map((column) => <th key={column}>{column}</th>)}</tr></thead>
                <tbody>
                  {rows.map((row, index) => (
                    <tr key={index}>{columns.map((column) => <td key={column}>{Array.isArray(row[column]) ? JSON.stringify(row[column]) : String(row[column] ?? '')}</td>)}</tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      )}
    </section>
  );
}

function PreviewPanel(props: { path: string; onPathChange: (path: string) => void }) {
  const [path, setPath] = React.useState(props.path);
  const [index, setIndex] = React.useState(0);
  const [summary, setSummary] = React.useState<UofsSummary | null>(null);
  const [summaryError, setSummaryError] = React.useState('');
  const [browserOpen, setBrowserOpen] = React.useState(false);
  const objectCount = summary?.objects ?? 0;
  const maxIndex = Math.max(0, objectCount - 1);
  React.useEffect(() => {
    setPath(props.path);
    setIndex(0);
  }, [props.path]);
  React.useEffect(() => {
    if (!path) {
      setSummary(null);
      setSummaryError('');
      return () => {};
    }
    let cancelled = false;
    setSummary(null);
    setSummaryError('');
    api.uofsSummary(path)
      .then((next) => {
        if (cancelled) return;
        setSummary(next);
        setIndex((current) => Math.min(current, Math.max(0, next.objects - 1)));
      })
      .catch((err) => {
        if (!cancelled) setSummaryError(err instanceof Error ? err.message : String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [path]);
  return (
    <section className="panel preview-panel">
      <div className="panel-header">
        <div>
          <h2>Visualize</h2>
          <p>Interactive 3D UofS object viewer.</p>
        </div>
      </div>
      <label className="field">
        <span>UofS path</span>
        <div className="path-control">
          <input value={path} onChange={(event) => { setPath(event.currentTarget.value); props.onPathChange(event.currentTarget.value); setIndex(0); }} />
          <button onClick={() => setBrowserOpen((open) => !open)}><Folder size={16} /> Browse data/</button>
        </div>
      </label>
      {browserOpen && <UofsBrowser currentPath={path} onClose={() => setBrowserOpen(false)} onChoose={(nextPath) => { setPath(nextPath); props.onPathChange(nextPath); setIndex(0); setBrowserOpen(false); }} />}
      {summaryError && <div className="errors">{summaryError}</div>}
      {!path ? (
        <EmptyState title="No file selected" detail="Please choose a SAP/2SAP file to visualize." />
      ) : (
        <>
          <div className="preview-controls">
            <button disabled={index <= 0 || objectCount === 0} onClick={() => setIndex(Math.max(0, index - 1))}>Previous</button>
            <div className="index-selector">
              <span>Object</span>
              <input
                type="number"
                value={objectCount > 0 ? index + 1 : ''}
                onChange={(event) => {
                  const val = event.currentTarget.value === '' ? 1 : Number(event.currentTarget.value);
                  setIndex(Math.max(0, Math.min(objectCount - 1, val - 1)));
                }}
              />
              <span>of {objectCount || 0}</span>
            </div>
            <button disabled={objectCount === 0 || index >= maxIndex} onClick={() => setIndex(Math.min(maxIndex, index + 1))}>Next</button>
          </div>
          <PolygonPreview path={path} index={index} />
        </>
      )}
    </section>
  );
}

function Toast(props: { message: string; onDismiss: () => void }) {
  React.useEffect(() => {
    if (!props.message) return undefined;
    const timer = window.setTimeout(props.onDismiss, 3500);
    return () => window.clearTimeout(timer);
  }, [props.message]);
  if (!props.message) return null;
  return <div className="toast"><span>{props.message}</span><button onClick={props.onDismiss}><X size={16} /></button></div>;
}

function App() {
  const [view, setView] = React.useState<'launcher' | 'jobs' | 'results' | 'visualize' | 'analysis'>('launcher');
  const [selectedJob, setSelectedJob] = React.useState<Job | undefined>();
  const [uofsPath, setUofsPath] = React.useState('');
  const [analysisTarget, setAnalysisTarget] = React.useState<{ action?: AnalysisAction }>({});
  const [toast, setToast] = React.useState('');
  function showToast(message: string) {
    setToast(message);
  }
  return (
    <div className="app-shell">
      <aside>
        <h1>SAP Workbench</h1>
        <button className={view === 'launcher' ? 'active' : ''} onClick={() => setView('launcher')}><Play size={17} /> Run Launcher</button>
        <button className={view === 'jobs' ? 'active' : ''} onClick={() => setView('jobs')}><Activity size={17} /> Jobs</button>
        <button className={view === 'results' ? 'active' : ''} onClick={() => setView('results')}><Database size={17} /> Data</button>
        <button className={view === 'visualize' ? 'active' : ''} onClick={() => setView('visualize')}><BarChart3 size={17} /> Visualize</button>
        <button className={view === 'analysis' ? 'active' : ''} onClick={() => { setView('analysis'); setAnalysisTarget({ action: 'summary' }); }}><ClipboardCheck size={17} /> Analysis</button>
      </aside>
      <main>
        {view === 'launcher' && <Launcher onToast={showToast} onJob={(job) => { setSelectedJob(job); setView('jobs'); }} />}
        {view === 'jobs' && <JobsPanel selectedJob={selectedJob} onSelect={setSelectedJob} onToast={showToast} onVisualize={(path) => { setUofsPath(path); setView('visualize'); }} onAnalyze={(path, action) => { setUofsPath(path); setAnalysisTarget({ action }); setView('analysis'); }} />}
        {view === 'results' && <DataPanel onToast={showToast} onVisualize={(path) => { setUofsPath(path); setView('visualize'); }} onAnalyze={(path, action) => { setUofsPath(path); setAnalysisTarget({ action }); setView('analysis'); }} />}
        {view === 'visualize' && <PreviewPanel path={uofsPath} onPathChange={setUofsPath} />}
        {view === 'analysis' && <AnalysisPanel path={uofsPath} action={analysisTarget.action} onPathChange={setUofsPath} onToast={showToast} />}
      </main>
      <Toast message={toast} onDismiss={() => setToast('')} />
    </div>
  );
}

createRoot(document.getElementById('root')!).render(<App />);
