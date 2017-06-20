"""
Get Data
"""
import iris
iris.FUTURE.netcdf_promote = True

def get_cube(inpath):
    print('In IO.get_cube')
    return iris.load(inpath)

def writeNetcdf(cb,outpath):
    print('In IO.writeNetcdf')
    return     
