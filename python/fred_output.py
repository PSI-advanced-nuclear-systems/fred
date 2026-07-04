"""
fred_output — shared HDF5 output utilities for the FRED 2.0 platform.

Provides:
  write_results_h5()   — dump solver results to an HDF5 file after run()
  read_results_h5()    — load results.h5 back as a dict of numpy arrays

Supports both FRED-ROD and FRED-OX solvers.

Units (FRED defaults, stored as-is):
  Temperature      : K
  Stress           : MPa
  Radii / lengths  : m
  Gap width        : m
  Burnup           : MWd/kgU
  Gas pressure     : MPa
  Fission gas      : mol
  Time             : s

HDF5 layout:
  /metadata          — attributes: app, nz, nf, nc, n_steps, all_steps
  /time              — [n_steps]                    seconds
  /thermal/T         — [n_steps, nz, nf+nc]         K
  /thermal/peak_T_fuel — [n_steps]                  K (peak fuel centreline)
  /thermal/gap_width   — [n_steps]                  m (axial-avg)
  /mechanical/       — FRED-ROD only
      pfc            — [n_steps]  MPa  pellet-cladding contact pressure
      rfo            — [n_steps]  m    fuel outer radius (axial avg)
      rci            — [n_steps]  m    cladding inner radius (axial avg)
      rco            — [n_steps]  m    cladding outer radius (axial avg)
      sigh_outer     — [n_steps]  MPa  outer clad hoop stress (axial avg)
  /burnup/           — FRED-OX only
      gpres          — [n_steps]  MPa  internal gas pressure
      fggen          — [n_steps]  mol  total fission gas generated
      fgrel          — [n_steps]  mol  total fission gas released
      bup            — [n_steps]  MWd/kgU  average burnup
  /restart/
      y              — [n_steps, neq]  IDA y-vector
      yp             — [n_steps, neq]  IDA yp-vector (time derivatives)
"""

import numpy as np

try:
    import h5py
except ImportError:
    raise ImportError(
        "h5py is required for HDF5 output. "
        "Install with: conda install h5py  or  pip install h5py"
    )


def write_results_h5(solver, nz, nf, nc, filename='results.h5', app='fred-rod',
                     all_steps=False):
    """
    Write solver results to an HDF5 file.

    Parameters
    ----------
    solver    : C++ FredRodSolver or FredOxSolver object (after run())
    nz, nf, nc: geometry parameters (number of axial/fuel/clad nodes)
    filename  : output filename (default 'results.h5')
    app       : 'fred-rod' or 'fred-ox'
    all_steps : whether all-steps mode was used (stored as metadata)
    """
    times = np.array(solver.time_points())
    n_steps = len(times)

    with h5py.File(filename, 'w') as f:

        # -- metadata --
        meta = f.create_group('metadata')
        meta.attrs['app']      = app
        meta.attrs['nz']       = nz
        meta.attrs['nf']       = nf
        meta.attrs['nc']       = nc
        meta.attrs['n_steps']  = n_steps
        meta.attrs['all_steps'] = int(all_steps)

        # -- time --
        f.create_dataset('time', data=times, compression='gzip')
        f['time'].attrs['units'] = 's'

        # -- thermal --
        th = f.create_group('thermal')
        T_raw = np.array(solver.temperatures())  # [n_steps, nz*(nf+nc)]
        T = T_raw.reshape(n_steps, nz, nf + nc)
        th.create_dataset('T', data=T, compression='gzip')
        th['T'].attrs['units'] = 'K'
        th['T'].attrs['shape'] = 'n_steps x nz x (nf+nc)'
        th['T'].attrs['description'] = (
            'Temperature at each node: [:, j, 0:nf] = fuel nodes (inner→outer), '
            '[:, j, nf:nf+nc] = cladding nodes (inner→outer) for axial layer j'
        )

        peak_T = np.array(solver.peak_fuel_temperature())
        th.create_dataset('peak_T_fuel', data=peak_T, compression='gzip')
        th['peak_T_fuel'].attrs['units'] = 'K'

        gap_w = np.array(solver.gap_width())
        th.create_dataset('gap_width', data=gap_w, compression='gzip')
        th['gap_width'].attrs['units'] = 'm'
        th['gap_width'].attrs['description'] = 'Axial-average radial gap width'

        # -- mechanical (FRED-ROD only) --
        if app == 'fred-rod':
            mech = f.create_group('mechanical')
            mech.create_dataset('pfc',
                data=np.array(solver.contact_pressure()), compression='gzip')
            mech['pfc'].attrs['units'] = 'MPa'
            mech['pfc'].attrs['description'] = 'Pellet-cladding contact pressure (axial avg)'

            mech.create_dataset('rfo',
                data=np.array(solver.fuel_outer_radius()), compression='gzip')
            mech['rfo'].attrs['units'] = 'm'

            mech.create_dataset('rci',
                data=np.array(solver.clad_inner_radius()), compression='gzip')
            mech['rci'].attrs['units'] = 'm'

            mech.create_dataset('rco',
                data=np.array(solver.clad_outer_radius()), compression='gzip')
            mech['rco'].attrs['units'] = 'm'

            mech.create_dataset('sigh_outer',
                data=np.array(solver.clad_outer_hoop_stress()), compression='gzip')
            mech['sigh_outer'].attrs['units'] = 'MPa'
            mech['sigh_outer'].attrs['description'] = 'Outer cladding hoop stress (axial avg)'

        # -- burnup / fission gas (FRED-OX only) --
        if app == 'fred-ox':
            bup_grp = f.create_group('burnup')
            bup_grp.create_dataset('gpres',
                data=np.array(solver.gas_pressure()), compression='gzip')
            bup_grp['gpres'].attrs['units'] = 'MPa'

            bup_grp.create_dataset('fggen',
                data=np.array(solver.fg_generated()), compression='gzip')
            bup_grp['fggen'].attrs['units'] = 'mol'
            bup_grp['fggen'].attrs['description'] = 'Total fission gas generated'

            bup_grp.create_dataset('fgrel',
                data=np.array(solver.fg_released()), compression='gzip')
            bup_grp['fgrel'].attrs['units'] = 'mol'
            bup_grp['fgrel'].attrs['description'] = 'Total fission gas released'

            bup_grp.create_dataset('bup',
                data=np.array(solver.burnup()), compression='gzip')
            bup_grp['bup'].attrs['units'] = 'MWd/kgU'
            bup_grp['bup'].attrs['description'] = 'Axial-average burnup'

        # -- restart (IDA y/yp vectors) --
        y_arr  = np.array(solver.y_out())   # [n_steps, neq]
        yp_arr = np.array(solver.yp_out())  # [n_steps, neq]
        if y_arr.size > 0:
            rst = f.create_group('restart')
            rst.create_dataset('y',  data=y_arr,  compression='gzip')
            rst.create_dataset('yp', data=yp_arr, compression='gzip')
            rst['y'].attrs['description']  = 'IDA state vector (y) at each output step'
            rst['yp'].attrs['description'] = 'IDA derivative vector (yp) at each output step'
            rst.attrs['note'] = (
                'Restart: reconstruct solver with same geometry/materials, '
                'then load y[-1] and yp[-1] as initial conditions.'
            )

    return filename


def read_results_h5(filename='results.h5'):
    """
    Read an HDF5 results file written by write_results_h5().

    Returns a dict with numpy arrays for all stored quantities.
    Keys depend on the app ('fred-rod' or 'fred-ox').

    Always present:
      time, T, peak_T_fuel, gap_width
      metadata: app, nz, nf, nc, n_steps, all_steps

    FRED-ROD only:
      pfc, rfo, rci, rco, sigh_outer

    FRED-OX only:
      gpres, fggen, fgrel, bup

    Restart (when IDA was active):
      restart_y, restart_yp
    """
    result = {}

    with h5py.File(filename, 'r') as f:
        # metadata
        meta = f['metadata']
        result['app']       = meta.attrs.get('app', 'fred-rod')
        result['nz']        = int(meta.attrs['nz'])
        result['nf']        = int(meta.attrs['nf'])
        result['nc']        = int(meta.attrs['nc'])
        result['all_steps'] = bool(meta.attrs.get('all_steps', 0))

        # time — use actual dataset length so a partial (crash-truncated) file still loads
        result['time'] = f['time'][:]
        result['n_steps'] = len(result['time'])

        # thermal
        # The C++ writer (FredSolverBase::openH5File/flushH5Step) stores T as
        # a flat 2-D dataset [n_steps, nz*(nf+nc)] — reshape to the documented
        # [n_steps, nz, nf+nc] layout expected by callers.
        T_raw = f['thermal/T'][:]
        if T_raw.ndim == 2:
            T_raw = T_raw.reshape(result['n_steps'], result['nz'], result['nf'] + result['nc'])
        result['T']           = T_raw
        result['peak_T_fuel'] = f['thermal/peak_T_fuel'][:]
        result['gap_width']   = f['thermal/gap_width'][:]

        # mechanical (fred-rod)
        if 'mechanical' in f:
            result['pfc']        = f['mechanical/pfc'][:]
            result['rfo']        = f['mechanical/rfo'][:]
            result['rci']        = f['mechanical/rci'][:]
            result['rco']        = f['mechanical/rco'][:]
            result['sigh_outer'] = f['mechanical/sigh_outer'][:]

        # burnup (fred-ox and fred-m-na)
        if 'burnup' in f:
            bup_grp = f['burnup']
            result['gpres'] = bup_grp['gpres'][:]
            result['fggen'] = bup_grp['fggen'][:]
            result['fgrel'] = bup_grp['fgrel'][:]
            result['bup']   = bup_grp['bup'][:]
            # fred-m-na extras
            for key in ('xwast', 'swtot', 'swopen', 'burst', 'melt'):
                if key in bup_grp:
                    result[key] = bup_grp[key][:]

        # restart
        if 'restart' in f:
            result['restart_y']  = f['restart/y'][:]
            result['restart_yp'] = f['restart/yp'][:]

    return result
