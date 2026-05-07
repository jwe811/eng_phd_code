import React from 'react';
import { createRoot } from 'react-dom/client';
import { Activity, BarChart3, Box, Database, FileText, Folder, Play, RefreshCw, StopCircle, X } from 'lucide-react';
import { api, commandPreview, DataBrowserListing, Job, JobRequest, ResultFile, UofsSummary, validateRequest } from './lib/api';
import { PolygonPreview } from './components/PolygonPreview';
import './styles.css';

const modes = ['Standard SAP', 'Hamiltonian SAP', '2SAP', '2SAP-Hamiltonian'];

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

function Launcher(props: { onJob: (job: Job) => void }) {
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
        <Play size={20} />
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

function JobsPanel(props: { selectedJob?: Job; onSelect: (job: Job) => void }) {
  const [jobs, setJobs] = React.useState<Job[]>([]);
  const [detail, setDetail] = React.useState<Job | undefined>(props.selectedJob);

  async function refresh() {
    const items = await api.jobs();
    setJobs(items);
    const selected = props.selectedJob ?? items[0];
    if (selected) setDetail(await api.job(selected.id));
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
  }

  return (
    <section className="panel">
      <div className="panel-header">
        <div>
          <h2>Jobs</h2>
          <p>SQLite-backed local history and logs.</p>
        </div>
        <button className="icon" onClick={refresh} aria-label="Refresh jobs"><RefreshCw size={18} /></button>
      </div>
      <div className="job-list">
        {jobs.map((job) => (
          <button key={job.id} className={detail?.id === job.id ? 'row selected' : 'row'} onClick={async () => { props.onSelect(job); setDetail(await api.job(job.id)); }}>
            <span>{job.tool}</span><strong>{job.status}</strong><small>{job.id}</small>
          </button>
        ))}
      </div>
      {detail && (
        <div className="detail">
          <div className="detail-actions">
            <code>{detail.command.join(' ')}</code>
            {detail.status === 'running' && <button className="danger" onClick={cancel}><StopCircle size={16} /> Cancel</button>}
          </div>
          <pre>{detail.log_tail || 'No logs yet.'}</pre>
        </div>
      )}
    </section>
  );
}

function ResultsPanel(props: { onPreview: (path: string) => void }) {
  const [results, setResults] = React.useState<ResultFile[]>([]);
  React.useEffect(() => {
    api.results().then(setResults).catch(console.error);
  }, []);
  const grouped = results.reduce<Record<string, ResultFile[]>>((acc, item) => {
    acc[item.group] = [...(acc[item.group] ?? []), item];
    return acc;
  }, {});
  return (
    <section className="panel">
      <div className="panel-header">
        <div>
          <h2>Results</h2>
          <p>Repo-local data files and metadata sidecars.</p>
        </div>
        <Database size={20} />
      </div>
      {Object.entries(grouped).map(([group, items]) => (
        <div key={group} className="result-group">
          <h3>{group}</h3>
          {items.slice(0, 80).map((item) => (
            <div className="result-row" key={item.path}>
              <div><strong>{item.path}</strong><small>{item.kind} · {item.size_bytes} bytes {item.meta_path ? '· meta' : ''}</small></div>
              {item.kind === 'uofs' && <button onClick={() => props.onPreview(item.path)}>Preview</button>}
            </div>
          ))}
        </div>
      ))}
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

function PreviewPanel(props: { path?: string }) {
  const [path, setPath] = React.useState(props.path ?? 'data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt');
  const [index, setIndex] = React.useState(0);
  const [summary, setSummary] = React.useState<UofsSummary | null>(null);
  const [summaryError, setSummaryError] = React.useState('');
  const [browserOpen, setBrowserOpen] = React.useState(false);
  const objectCount = summary?.objects ?? 0;
  const maxIndex = Math.max(0, objectCount - 1);
  React.useEffect(() => {
    if (props.path) {
      setPath(props.path);
      setIndex(0);
    }
  }, [props.path]);
  React.useEffect(() => {
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
          <h2>Preview</h2>
          <p>Interactive 3D UofS object viewer.</p>
        </div>
        <Box size={20} />
      </div>
      <label className="field">
        <span>UofS path</span>
        <div className="path-control">
          <input value={path} onChange={(event) => { setPath(event.currentTarget.value); setIndex(0); }} />
          <button onClick={() => setBrowserOpen((open) => !open)}><Folder size={16} /> Browse data/</button>
        </div>
      </label>
      {browserOpen && <UofsBrowser currentPath={path} onClose={() => setBrowserOpen(false)} onChoose={(nextPath) => { setPath(nextPath); setIndex(0); setBrowserOpen(false); }} />}
      {summaryError && <div className="errors">{summaryError}</div>}
      <div className="preview-controls">
        <button disabled={index <= 0 || objectCount === 0} onClick={() => setIndex(Math.max(0, index - 1))}>Previous</button>
        <span>{objectCount > 0 ? `Object ${index + 1} of ${objectCount}` : 'No objects loaded'}</span>
        <button disabled={objectCount === 0 || index >= maxIndex} onClick={() => setIndex(Math.min(maxIndex, index + 1))}>Next</button>
      </div>
      <PolygonPreview path={path} index={index} />
    </section>
  );
}

function App() {
  const [view, setView] = React.useState<'launcher' | 'jobs' | 'results' | 'preview'>('launcher');
  const [selectedJob, setSelectedJob] = React.useState<Job | undefined>();
  const [previewPath, setPreviewPath] = React.useState<string | undefined>();
  return (
    <div className="app-shell">
      <aside>
        <h1>SAP Workbench</h1>
        <button className={view === 'launcher' ? 'active' : ''} onClick={() => setView('launcher')}><Play size={17} /> Run Launcher</button>
        <button className={view === 'jobs' ? 'active' : ''} onClick={() => setView('jobs')}><Activity size={17} /> Jobs</button>
        <button className={view === 'results' ? 'active' : ''} onClick={() => setView('results')}><Database size={17} /> Results</button>
        <button className={view === 'preview' ? 'active' : ''} onClick={() => setView('preview')}><BarChart3 size={17} /> Preview</button>
      </aside>
      <main>
        {view === 'launcher' && <Launcher onJob={(job) => { setSelectedJob(job); setView('jobs'); }} />}
        {view === 'jobs' && <JobsPanel selectedJob={selectedJob} onSelect={setSelectedJob} />}
        {view === 'results' && <ResultsPanel onPreview={(path) => { setPreviewPath(path); setView('preview'); }} />}
        {view === 'preview' && <PreviewPanel path={previewPath} />}
      </main>
    </div>
  );
}

createRoot(document.getElementById('root')!).render(<App />);
