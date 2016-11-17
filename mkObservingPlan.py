import astroplan
from astroplan import Observer, FixedTarget, ObservingBlock
from astroplan.constraints import AtNightConstraint, AirmassConstraint,\
                    TimeConstraint
from astroplan.scheduling import Transitioner, SequentialScheduler, Schedule,\
                    PriorityScheduler
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astropy import units as u
import pandas as pd
import numpy as np
from astLib import astCoords


def mst2utc(time):
    ''' Converts from MST to UTC '''
    dt = TimeDelta(7*3600, format='sec')
    return time + dt

def utc2mts(time):
    ''' Converts from UTC to MST '''
    dt = TimeDelta(7*3600, format='sec')
    return time - dt


# load the data from the CSV of the spreadsheet
data = pd.read_csv('./PSZ2_unconfirmed_catalog.csv')

# figure out which we need to observe
mask = pd.isnull(data['IR run1'])
data = data[mask]

mask = data['SNR'] > 5.
data = data[mask]

# make telescope location
kpno = Observer.at_site('KPNO', timezone='US/Mountain')

# make target coordinates
coords = [SkyCoord(ra, dec, unit='deg', frame='icrs') for ra, dec in
          zip(data['RA'], data['DEC'])]

# make all the targets
targets = [FixedTarget(name=name, coord=coor) for name, coor in
           zip(data['Name'], coords)]

# make the observing time - Local time.
start_time_mst = Time('2016-11-14 1:30')
end_time_mst = Time('2016-11-14 07:00')
night_time = end_time_mst - start_time_mst
observable_time_mst = start_time_mst + night_time * np.linspace(0,1,75)

# convert everything to UTC
start_time_utc = mst2utc(start_time_mst)
end_time_utc = mst2utc(end_time_mst)
observable_time_utc = mst2utc(observable_time_mst)

# now we figure out what is observable at all with some simple constraints
constraint = [AtNightConstraint.twilight_civil(),
              AirmassConstraint(max=2.7, boolean_constraint=False)]

#observable = astroplan.is_observable(constraint, kpno, targets,
#                              times=observable_time_utc)

#data['Observable'] = observable
#data.to_csv('updated.csv')

# now we make the exposure times
exp_time_single = 12 * u.second
readout_time = 1.7 * u.second
n_coadd = 5
n_dither = 12
exp_time_tot = n_coadd * n_dither * exp_time_single + readout_time

no_optical_data = pd.isnull(data['Optical run1'])

# RA, DEC for the CCD positions -- in arcseconds
ccd_shift = [[427, 427], [0, -850], [-850, 0], [0, 850]]

# observing blocks
blocks = []


array = data['SNR'].values
temp = array.argsort()[::-1] # sort high to low
ranks = np.empty(len(array), int)
ranks[temp] = np.arange(len(array)) +1

for index, t in enumerate(targets):
#    if no_optical_data.values[index]:
#        priority = 1
#    else:
#        priority = 2

    b = ObservingBlock(t,
                        4*1230*u.second,
                        priority = ranks[index],
                        constraints=constraint,
                       configuration={'SNR':'%.4f' % data['SNR'].values[index]})
    blocks.append(b)

# now we make a transitioner
slew_rate = .5 * u.deg/u.second
setup_time = {'default': 160 * u.second}
transitioner = Transitioner(slew_rate)

'''
# now we do the scheduling
seq_scheduler = SequentialScheduler(constraints=constraint, observer=kpno,
                                    transitioner=transitioner)
sequential_schedule = Schedule(start_time_utc, end_time_utc)

try:
    seq_scheduler(blocks, sequential_schedule)
except ValueError:
    print(sequential_schedule.to_table())
'''

# Initialize the priority scheduler with the constraints and transitioner
prior_scheduler = PriorityScheduler(constraints=constraint,
                                     observer=kpno,
                                     transitioner=transitioner)
# Initialize a Schedule object, to contain the new schedule
priority_schedule = Schedule(start_time_utc, end_time_utc)

try:
    prior_scheduler(blocks, priority_schedule)
except ValueError:
    print(priority_schedule.to_table())

