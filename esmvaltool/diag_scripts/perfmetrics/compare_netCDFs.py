



import xarray as xr
import numpy as np
import pandas as pd


# Load the NetCDF files
file_diagnostics_path = '/home/b/b309265/ta850-global_to_sm-global_RMSD.nc'
file_performance_path = '/home/b/b309265/performance_metrics.nc'

# Load datasets
ds_diagnostics = xr.open_dataset(file_diagnostics_path)
ds_performance = xr.open_dataset(file_performance_path)

# Strip the "CMIP5_" prefix from model names in the performance dataset
ds_performance['alias'] = [model.replace("CMIP5_", "").strip() for model in ds_performance['alias'].values]

# Convert all model names and strip whitespace for consistency
models_diagnostics = ds_diagnostics['models'].astype(str).str.strip()
models_performance = ds_performance['alias'].astype(str).str.strip()

# Display unique variables from each dataset
variables_diagnostics = set(ds_diagnostics['diagnostics'].astype(str).values)
variables_performance = set(ds_performance['variable'].astype(str).values)

print("Variables in Diagnostics Dataset:", variables_diagnostics)
print("Variables in Performance Dataset:", variables_performance)

# Standardize variable names by converting to lowercase and stripping whitespace
variables_diagnostics = {var.lower().strip() for var in variables_diagnostics}
variables_performance = {var.lower().strip() for var in variables_performance}

# Find common models and standardized variables
common_models = np.intersect1d(models_diagnostics, models_performance)
common_variables = variables_diagnostics.intersection(variables_performance)

print("Common Models:", common_models)
print("Common Variables after Standardization:", common_variables)

# Initialize a dictionary to collect comparison results
comparison_results = {
    "Model": [],
    "Variable": [],
    "Grade_Value": [],
    "Var_Value": [],
    "Difference": []
}

# Perform comparison for each common model-variable pair
for model in common_models:
    for variable in common_variables:
        try:
            # Select values for each matched model-variable pair
            # Align 'diagnostics' and 'variable' names by standardizing to lowercase
            grade_val = ds_diagnostics.grade.sel(models=model, diagnostics=variable).values
            var_val = ds_performance.var.sel(alias=model, variable=variable).values

            # Calculate the absolute difference
            diff = np.abs(grade_val - var_val).flatten()
            
            # Store results
            comparison_results["Model"].append(model)
            comparison_results["Variable"].append(variable)
            comparison_results["Grade_Value"].append(grade_val)
            comparison_results["Var_Value"].append(var_val)
            comparison_results["Difference"].append(diff)
        
        except KeyError:
            # Skip if the pair is missing in one dataset
            print(f"Skipped {model}-{variable} due to mismatch in one of the datasets")

# Convert to DataFrame and check results
comparison_df = pd.DataFrame(comparison_results)
print("Comparison Results:\n", comparison_df)

# Optionally, save results to CSV
# comparison_df.to_csv("comparison_results.csv", index=False)
