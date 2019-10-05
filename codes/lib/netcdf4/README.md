# Purpose
A sanitizing layer for the netCDF4 library. Adds a number of convenience methods and aims for a cleaner
user interface. All classes avaliable are children of their associated netCDF4 counterparts.
The module is therefore usable as a drop-in replacement for the classes Dataset, Group, Dimension and
Variable of the netCDF4 module.

#Usage

The class NcDataset offers a context manager. All reading from and writing to a file can be moved
into the context manager block. Opening and closing will be done automatically:
```python
from netcdf4 import NcDataset

# open a file in write mode, create a dimension and variable
with NcDataset("test.nc", "w") as nc:
	nc.createDimension("time", None)
	var = nc.createVariable("time", "i32", ("time",))
	var[:] = range(100)
# the file will be closed by now

# open a file to read data
with NcDataset("test.nc", "r") as nc:
	data = nc.variables["time"][:]
```

NcDataset, NcGroup and NcVariables all offer copy* methods to easily transfer content from one
file to another:
```python
# copy content from one file to another
with NcDataset("infile.nc", "r") as ncin:
	with NcDataset("outfile.nc", "w") as ncout:
		# copy the entire content of ncin into ncout
		ncout.copyDataset(ncin)
		# copy the content with out specific dimension/variable/group
		ncout.copyDataset(ncin, skipdims=("x","y"), skipvars=("x","y"), skipgroups=("group1"))
		# copy only the variable x from ncin to ncout 
		ncout.copyVariable(ncin["x"])
		# copy only the declaration of variable x from ncin to ncout, without copying the data
		ncout.copyVariable(ncin["x"], data=False)
		# copy all variables from ncin to ncout
		ncout.copyVariables(ncin.variables)
		# copy all but the variables x and y from ncin to ncout
		ncout.copyVariables(ncin.variables, skip=("x","y"))
		# same as above but without copying the data
		ncout.copyVariables(ncin.variables, skip=("x","y"), data=False)
		# similar routines exist from groups dimensions and attributes
        #
		# copy whole dataset but one variable that is changed in the code
        ncout.copyDataset(ncin, skipvars="var_name")
        var = ncout.copyVariable(ncin["var_name"], data=False)
        var[:] = some_local_variable
```

netcdf4 offers some filtering functionality. NcDataset/NcGroup classes offer filter* methods:
```python
with NcDataset("test.nc", "r") as nc:
	# get all dimensions with a length of 10
	nc.filterDimensions(10)
	# get all variables with exactly 3 dimensions
	nc.filterVariables(ndim=3)
	# get all variables depending on the dimension time and x
	nc.filterVariables(dims=("time", "x"))
```

The module offers some API consolidation with regard to attributes: 
```python
with NcDataset("test.nc", "w") as nc:
	# create or overwrite the attribute unit with value "mm/d"
	nc.setAttribute("unit", "mm/d")
	# create or overwrite multiple attributes
	nc.setAttributes({"unit": "mm/d", "long_name": "precipitation"})
	# retrieve all attributes
	nc.getAttributes()
```

Some random query methods:
```python
with NcDataset("test.nc", "r") as nc:
	# get all dates as datetime.datetime objects
	nc.getDates()
	# get the fill value
	nc.getFillValue()
```
