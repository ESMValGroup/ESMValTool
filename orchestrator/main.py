from mocking import io, backend
inpath="test/testdata/DummyM3.nc"
outpath="/tmp/testEVT/"

cb = io.get_cube(inpath)
cb = backend.timeExtract(cb)
cb = backend.reformat(cb)
cb = backend.extractLevel(cb)
cb = backend.regrid(cb)
cb = backend.mask(cb)
cb = backend.gridPointOperation(cb)
cb = backend.timeOperations(cb)
cb = backend.multiModelstatistics(cb)
io.writeNetcdf(cb,outpath)
