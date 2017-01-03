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
    dt = TimeDelta(7 * 3600, format='sec')
    return time + dt


def utc2mts(time):
    ''' Converts from UTC to MST '''
    dt = TimeDelta(7 * 3600, format='sec')
    return time - dt

# load the data from the CSV of the spreadsheet
data = pd.read_csv('./catalogs/PSZ2_unconfirmed_catalog.csv')

# figure out which we need to observe
mask = pd.isnull(data['Optical run1'])
data = data[mask]

mask = data['SNR'] > 5.
data = data[mask]

# make telescope location
kpno = Observer.at_site('KPNO', timezone='US/Mountain')

# make target coordinates
coords = [SkyCoord(ra, dec, unit='deg', frame='icrs')
          for ra, dec in zip(data['RA'], data['DEC'])]

# make all the targets
targets = [FixedTarget(name=name, coord=coor)
           for name, coor in zip(data['Name'], coords)]

# make the observing time - Local time.
start_time_mst = Time('2016-11-23 18:30')
end_time_mst = Time('2016-11-26 07:00')
night_time = end_time_mst - start_time_mst
observable_time_mst = start_time_mst + night_time * np.linspace(0, 1, 75)

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
exp_time_single_gr = 90 * u.second
exp_time_single_iz = 275 * u.second
readout_time = 30 * u.second
settle_time = 10 * u.second
n_coadd = 1
n_dither = 4
exp_time_tot_gr = n_coadd * n_dither * (
    exp_time_single_gr + readout_time + settle_time)

exp_time_tot_iz = n_coadd * n_dither * (
    exp_time_single_iz + readout_time + settle_time)

array = data['SNR'].values
temp = array.argsort()[::-1]  # sort high to low
ranks = np.empty(len(array), int)
ranks[temp] = np.arange(len(array)) + 1

# observing blocks
blocks = []
for index, t in enumerate(targets):
    b = ObservingBlock(t,
                       3600 * u.second,
                       priority=ranks[index],
                       constraints=constraint,
                       configuration={'SNR':
                                      '%.4f' % data['SNR'].values[index]})
    blocks.append(b)

# now we make a transitioner
slew_rate = .5 * u.deg / u.second
transitioner = Transitioner(slew_rate)

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
