import numpy as np
import LYRA_functions as lyra
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import time

ra_input = np.random.uniform(0.00001,359.111,int(5e3))
dec_input = np.random.uniform(-89.999,87.000,int(5e3))
df = pd.DataFrame(data={'ra': ra_input, 'dec': dec_input})
# id=np.random.uniform(0,1e4,int(1e4))

start = time.time()
names = lyra.LYRA_source_names(ra_input, dec_input)
end = time.time()
# c = SkyCoord(ra=ra_input*u.degree, dec=dec_input*u.degree)
# cc = c.to_string('hmsdms')

# df['names'] = df.apply( deg2sex(df['ra'], df['dec']) )
print(end-start)
print(type(names))