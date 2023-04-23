## <img src="/Users/suzhenbo/Mylibrary/Projects/jupyter_plot/disk/jl_version/AGNDisk/test/logo.png" alt="logo" style="zoom:15%;" />Accretion Disk Code

### General views

- [x] Classical SSD model, [Shakura & Sunyaev 1973](https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract) .

- [x] Line-driven wind disk model, [Laor & Davis 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.3024L/abstract).

### Basic usage

```julia
using Pkg
activate ./AGNDisk
#using AGNDisk
push!(LOAD_PATH, "./jl_version/AGNDisk/src/")
using AGNDisk
```

```julia
p = DiskPars(bhmass=mbh, mdot=mdot, aspin=aspin) # Initial of disk parameters.
t = thindisk_sed(p, wavelength) # Modeling SSD model.
w = winddisk_sed(p, wavelength) # Modeling line-driven disk model.

save_res(t, p; fi_np="./sim_lumin") # Save results.
save_res(w, p; fi_np="./sim_lumin") # Save results.
```

```julia
showdata(p; fi_np="./sim_lumin") # Plot of data
showdata(w; fi_np="./sim_lumin") # Plot of data
```

![thindisk](/Users/suzhenbo/Mylibrary/Projects/jupyter_plot/disk/jl_version/AGNDisk/test/thindisk.png)![winddisk](/Users/suzhenbo/Mylibrary/Projects/jupyter_plot/disk/jl_version/AGNDisk/test/winddisk.png)

### TEST

- For more info please see [Laor & Davis 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.3024L/abstract).

![laor](/Users/suzhenbo/Mylibrary/Projects/jupyter_plot/disk/jl_version/AGNDisk/test/laor.png)



