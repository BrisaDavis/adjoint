"""
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

""" 

import os
import numpy as np


#-----------------------------------------------
# Set these parameters for adjoint flagging....

# location of output from computing adjoint:
adjoint_output = os.path.abspath('../adjoint/_output')
print('Will flag using adjoint solution from  %s' % adjoint_output)

# Time period of interest:
t1 = 3.5*3600.
t2 = 11*3600.

# tolerance for adjoint flagging:
adjoint_flag_tolerance = 0.004
#-----------------------------------------------

t_shelf = 3.2*3600   # time approaching continental slope
t_harbor = 3.5*3600  # time approaching harbor

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun. 

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    from clawpack.clawutil import data 
    
    
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Adjoint specific data:
    #------------------------------------------------------------------
    rundata = setadjoint(rundata)
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 140.0          # xlower
    clawdata.upper[0] = 250.0          # xupper
    clawdata.lower[1] = 10.0          # ylower
    clawdata.upper[1] = 62.0          # yupper
    
    # Number of grid cells:
    clawdata.num_cells[0] = 110      # mx
    clawdata.num_cells[1] = 52      # my
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 4
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 1
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 22
        clawdata.tfinal = 11*3600.
        clawdata.output_t0 = False # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  list(np.linspace(3600,3600*9,9))
 
    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format = 'binary'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 1
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer']
    
    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper
                  
       
    # ---------------
    # gauges:
    # ---------------

    gauges = rundata.gaugedata.gauges
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    # Outside harbor:
    gauges.append([1, 235.536, 41.67, t_shelf, 1.e10])
        
    # Inside harbor:
    gauges.append([2, 235.80917,41.74111,t_harbor, 1.e10])

                  
    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

    

    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [5, 6, 6, 3, 30]
    amrdata.refinement_ratios_y = [5, 6, 6, 3, 30]
    amrdata.refinement_ratios_t = [5, 6, 6, 3, 4]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    # need 4 values, set in setadjoint


    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    # see setadjoint to set tolerance for adjoint flagging

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3       

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0      


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions 
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    regions.append([1, 1, 0., 1e9, 0, 360, -90, 90]) #whole world
    regions.append([1, 3, 0., 7*3600., 0, 360, -90, 90]) #whole world
    regions.append([1, 3, 7*3600.,10*3600., 170., 360, 18, 90])
    regions.append([1, 3, 10*3600.,1e9, 195., 360, -90, 90])
    regions.append([4, 4, 0., 1800, 175, 195, 50, 54]) #earthquake source AASZ04

    regions.append([3, 4, t_shelf, 1e9, 235, 238, 34, 43]) # between shelf and CC
    regions.append([4, 4, t_shelf, 1e9, 235, 236, 41, 42]) 
    regions.append([5, 5, t_shelf, 1e9, 235.5,235.83,41.6,41.8]) #only harbor 
    regions.append([5, 6, t_harbor, 1e9, 235.78,235.84,41.735,41.775]) #only harbor

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata

    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geo_data attribute"
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    tide_stage = 77.
    geo_data.sea_level = (tide_stage - 77.)/100.    #  m relative to MHW

    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 100.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.09
    refinement_data.deep_depth = 100.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    
    topo_path = os.path.join(scratch_dir, 'etopo1min170E124W40N61N.asc')
    topo_data.topofiles.append([3, 1, 1, 0., 1.e10, topo_path])
    
    topo_path = os.path.join(scratch_dir, 'etopo4min120E110W0N62N.asc')
    topo_data.topofiles.append([3, 1, 1, 0., 1.e10, topo_path])
    
    topo_path = os.path.join(scratch_dir, 'cc-1sec-c.asc')
    topo_data.topofiles.append([-3, 1, 1, 32000., 1.e10, topo_path])
    
    topo_path = os.path.join(scratch_dir, 'cc-1_3sec-c_pierless.asc')
    topo_data.topofiles.append([3, 1, 1, 32000., 1.e10, topo_path])


    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :
    #   [topotype, minlevel,maxlevel,fname]
    dtopo_path = os.path.join(scratch_dir, 'AASZ04v2.tt3')
    dtopo_data.dtopofiles.append([3,3,3,dtopo_path])
    dtopo_data.dt_max_dtopo = 0.2


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  0
    rundata.qinit_data.qinitfiles = []
    qinitfiles = rundata.qinit_data.qinitfiles 
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == fixedgrids.data values ==
    rundata.fixed_grid_data.fixedgrids = []
    fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files

    return rundata
    # end of function setgeo
    # ----------------------

#-------------------
def setadjoint(rundata):
#-------------------

    """
    Set parameters used for adjoint flagging.
    Also reads in all of the checkpointed Adjoint files.
    """
        
    import glob

    # Set these parameters at top of this file:
    # adjoint_flag_tolerance, t1, t2, adjoint_output
    # Then you don't need to modify this function...

    rundata.amrdata.flag2refine = True  # for adjoint flagging
    rundata.amrdata.flag2refine_tol = adjoint_flag_tolerance
    
    rundata.clawdata.num_aux = 4
    rundata.amrdata.aux_type = ['center', 'capacity', 'yleft', 'center']
    
    adjointdata = rundata.new_UserData(name='adjointdata',fname='adjoint.data')
    adjointdata.add_param('adjoint_output',adjoint_output,'adjoint_output')
    adjointdata.add_param('t1',t1,'t1, start time of interest')
    adjointdata.add_param('t2',t2,'t2, final time of interest')
    
    files = glob.glob(os.path.join(adjoint_output,"fort.b*"))
    files.sort()
    
    if (len(files) == 0):
        print("No binary files found for adjoint output!")
    
    adjointdata.add_param('numadjoints', len(files), 'Number of adjoint output files.')
    adjointdata.add_param('innerprod_index', 4, 'Index for innerproduct data in aux array.')
    
    counter = 1
    for fname in files:
        f = open(fname)
        time = f.readline().split()[-1]
        adjointdata.add_param('file' + str(counter), fname, 'Binary file' + str(counter))
        counter = counter + 1

    return rundata
    # end of function setadjoint
    # ----------------------

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
