
# coding: utf-8

# In[244]:

import iris
import os
import glob
import datetime
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np


# In[254]:

data_path = 'ftp.pa.op.dlr.de/pub/righi/BACKEND-DATA/ETHZ_CMIP5/historical/Amon'
ps_cube_list = iris.load(glob.glob(os.path.join(data_path, "ps/GFDL-ESM2G/r1i1p1/ps_Amon_GFDL-ESM2G_historical_r1i1p1_*")))
tro3_cube = iris.load_cube(os.path.join(data_path, "tro3/GFDL-ESM2G/r1i1p1/tro3_Amon_GFDL-ESM2G_historical_r1i1p1_200001-200212.nc"))

for cube in ps_cube_list:
    cube.attributes['history'] = None
    
equalise_attributes(ps_cube_list)
ps_cube = ps_cube_list.concatenate_cube()
print ps_cube.shape
print tro3_cube.shape


# In[246]:

def extract_time_period(cube, start_time, end_time):
    """
    Extract only the data between `start_time` and `end_time` from provided cube. 

    It uses the time point of the data, or if available the beginning of a time 
    period. This time point has to be at or after `start`, and 
    before or at `end`. 

    :param Cube cube: Iris Cube. 
    :param datetime.datetime start_time: 
    :param datetime.datetime end_time: 
    :returns: Cube containing only data between `start_time` and `end_time`. 
    :rtype: Iris Cube
    """
    if cube.coord('time').has_bounds():
        time_range_constraint = iris.Constraint(time = lambda cell: start_time <= cell.bound[0] and cell.bound[1] <= end_time)
    else:
        time_range_constraint = iris.Constraint(time = lambda cell: start_time <= cell.point <= end_time)

    with iris.FUTURE.context(cell_datetime_objects=True):
        return cube.extract(time_range_constraint) 


# start_time = datetime.datetime(2001, 1, 1)
# end_time = datetime.datetime(2003, 1, 1)
# 
# tro3_cube = extract_time_period(tro3_cube, start_time, end_time)

# In[248]:

print(ps_cube)
print(tro3_cube)


# In[250]:

l = tro3_cube.coord('air_pressure').points
l


# In[56]:

import matplotlib.pyplot as plt
import iris.quickplot as qplt

for i, p_lev in enumerate(tro3_cube.coord('air_pressure').points):
    print('Pressure level:' + str(p_lev))
    qplt.contourf(tro3_cube[0, i,:, :])

    # Add coastlines to the map
    plt.gca().coastlines()
    plt.show()


# In[342]:

def bounds_widths(values):
    """
    For each value in a 1D array, returns a 1D array with the bounds widths (the distances 
    between the bounds around each value).
    """
    assert iris.util.monotonic(values)
    
    values = np.append(values, values[-1])
    values = np.insert(values, 0, values[0])
    bounds_widths = []
    num_values = len(values)
    for i, val in enumerate(values):
        if i == 0 or i == num_values - 1:  # boundary values
            continue
        else:
            dist_to_lower_bound = values[i-1] - val
            dist_to_upper_bound = val - values[i+1]
            if i == 1:                     # second value
                bounds_width =       dist_to_lower_bound + 0.5 * dist_to_upper_bound
            elif i == num_values - 2:      # penultimate value
                bounds_width = 0.5 * dist_to_lower_bound +       dist_to_upper_bound
            else:
                bounds_width = 0.5 * dist_to_lower_bound + 0.5 * dist_to_upper_bound
        bounds_widths.append(bounds_width)
    return np.abs(bounds_widths)


values = np.array(range(10, 100, 10))
print values
print values.shape
print bounds_widths(values)
print bounds_widths(values).shape
print


# In[341]:

def pressure_layer_thickness(cube, top_limit=500):
    """ 
    """
    pressure_values = cube.coord('air_pressure').points
    #pressure_values = range(100, 10, -10)
    
    assert np.min(pressure_values) > top_limit
    
    #print pressure_values
    
    pressure_values = np.append(pressure_values, top_limit)
    print pressure_values
    
    p_diff = bounds_widths(pressure_values)
    p_diff = np.array(p_diff)
    #print p_diff.shape

    return p_diff
    
print(pressure_layer_thickness(tro3_cube, top_limit=500))


# In[ ]:

def pressure_layer_thickness_simple(cube, top_limit=5, bottom_limit=105):
    """ 
    """
    pressure_values = cube.coord('air_pressure').points
    #pressure_values = range(100, 10, -10)
    
    assert np.min(pressure_values) > top_limit
    assert np.min(pressure_values) < bottom_limit
    
    print pressure_values
    
    pressure_values = np.insert(pressure_values, 0, bottom_limit)
    pressure_values = np.append(pressure_values, top_limit)
    print pressure_values
    
    p_diff = bounds_widths(pressure_values)
    p_diff = np.array(p_diff)
    #print p_diff.shape

    return p_diff
    
print(pressure_layer_thickness(tro3_cube, top_limit=5, bottom_limit=105000))


# In[314]:

def pressure_layer_thickness_surface(tro3_cube, ps_cube, top_limit=100):
    pressure_4d = create_pressure_cube(tro3_cube, ps_cube)
    
    res = apply_bounds_widths(pressure_4d)
    
    return res    
    
    
def create_pressure_cube(tro3_cube, ps_cube):
    """
    Returns: A numpy array
    """
    #print tro3_cube.shape
    p = tro3_cube.coord('air_pressure').points

    p_cube = iris.util.broadcast_to_shape(p, tro3_cube.shape, [1])
    #print p_cube.shape

    ps_cube_bc = iris.util.broadcast_to_shape(ps_cube.data, tro3_cube.shape, [0, 2, 3])
    #print ps_cube_bc.shape

    pressure_4d = np.where((ps_cube_bc - p_cube) < 0, ps_cube_bc, p_cube)
    return pressure_4d


def apply_bounds_widths(cube):
    return np.apply_along_axis(bounds_widths, 1, cube)


# In[ ]:

def pressure_layer_thickness_surface(tro3_cube, ps_cube, top_limit=100):
    p_dim = 1 # TODO
    
    # TODO: Constrain to upperlimit=5mb
    
    pressure_4d = create_pressure_cube(tro3_cube, ps_cube)
    
    # TOD: Add top(5mb) and bottom(surface pressure) layer
    
    upper_contribution = 0.5 * (pressure_4d[:, 2:] - pressure_4d[:,1:-1])
    lower_contribution = 0.5 * (pressure_4d[:,1:-1] - pressure_4d[:, :-2])
    
    res = upper_contribution + lower_contribution
    
    res[:, 0]
    
    return res    


# In[313]:

pressure_4d = create_pressure_cube(tro3_cube, ps_cube)

# plot pressure slices
for i in range(pressure_4d.shape[1]):
    plt.contourf(pressure_4d[0, i,:, :])

    # Add coastlines to the map
    plt.show()


# In[322]:

# apply bounds_widths to pressure_4d
get_ipython().magic(u'time res = apply_bounds_widths(pressure_4d)')
print res.shape

# plot bounds_widths
for i in range(res.shape[1]):
    plt.contourf(res[0, i,:, :])

    # Add coastlines to the map
    plt.colorbar()
    plt.show()


# ![Ozone vs height](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Ozone_altitude_UV_graph.svg/300px-Ozone_altitude_UV_graph.svg.png)
# 
# ![Pressure vs Height](https://upload.wikimedia.org/wikipedia/commons/thumb/9/95/Pressure_air.svg/300px-Pressure_air.svg.png)

# In[293]:

# TODO: result units

g = 9.81  # m s^-2
mw_air = 29  # g/mol
mw_ozone = 48  # g/mol

def total_column_ozone(tro3_cube):
    pressure_layer_thickness = pressure_layer_thickness_simple(tro3_cube, top_limit=100, bottom_limit=102000)
    pressure_layer_thickness = pressure_layer_thickness[np.newaxis, :, np.newaxis, np.newaxis]
    return (tro3_cube * pressure_layer_thickness / g * mw_ozone / mw_air).collapsed('air_pressure', iris.analysis.SUM)
    #  mol_O3/mol_air       N m^-2               N kg^-1   g_O3/mol_O3   g_air/mol_air
    #                         m^-2                 kg      g_O3         g_air

ozone_col_cube = total_column_ozone(tro3_cube)
print ozone_col_cube

qplt.contourf(ozone_col_cube[0,:,:])
plt.gca().coastlines()
plt.title('Total column ozone')
plt.show()


# In[321]:

def total_column_ozone(tro3_cube, ps_cube):
    pressure_layer_thickness = res #pressure_layer_thickness_surface(tro3_cube, ps_cube, top_limit=100)
    return (tro3_cube * pressure_layer_thickness / g * mw_ozone / mw_air).collapsed('air_pressure', iris.analysis.SUM)
    #  mol_O3/mol_air       N m^-2               N kg^-1   g_O3/mol_O3   g_air/mol_air
    #                         m^-2                 kg      g_O3         g_air

ozone_col_cube = total_column_ozone(tro3_cube, ps_cube)
print ozone_col_cube

qplt.contourf(ozone_col_cube[0,:,:])
plt.gca().coastlines()
plt.title('Total column ozone')
plt.show()


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[319]:

a = np.array([0.0, 10.0, 20.0, 30.0])
b = np.array([1.0, 2.0, 3.0])
a[:, np.newaxis] + b


# In[ ]:



