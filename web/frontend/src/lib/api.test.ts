import { describe, expect, it } from 'vitest';
import { commandPreview, validateRequest } from './api';

describe('launcher validation', () => {
  it('rejects lattice products at ten or above', () => {
    const errors = validateRequest({
      tool: 'creator',
      L: 4,
      M: 1,
      mode: 0,
      span: 2,
      options: { export_eigenvectors: false, damping: false }
    });
    expect(errors.join(' ')).toContain('(L + 1)(M + 1)');
  });

  it('shows the command preview for creator runs', () => {
    const preview = commandPreview({
      tool: 'creator',
      L: 1,
      M: 1,
      mode: 0,
      span: 2,
      options: { export_eigenvectors: false, damping: false }
    });
    expect(preview).toBe('bin/creator_all -L 1 -M 1 -m 0 -s 2');
  });
});

