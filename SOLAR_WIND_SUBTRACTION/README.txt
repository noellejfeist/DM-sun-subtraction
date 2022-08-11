SOLAR WIND SUBTRACTION

Written by Noelle Feist
2022 Cornell University Astrophysics REU
Funded by NSF award NST/AST-1950324

Thanks to Shami Chatterjee, Thankful Cromartie, and Jim Cordes for helping me this summer.

--Files--
    solarWindPlots.py: The main script that creates all of the plots.
    solarWindTools.py: A helper script for solarWindPlots.py.
    parReader.py: A script that extracts data from a .par file and stores them in a custom Pulsar class.
    parToDmx.py: A script that uses parReader.py and exports a specified .par file as a .dmx file.

--solarWindPlots.py--
    Command Line: python3 solarWindPlots.py <input_dir>

    Usage: For each .par file in <input_dir>, the solar wind component of the dD (DMX) measurements will be removed. The sun-subtracted .dmx files are placed in a 
        subdirectory of <input_dir> called "dmx". Various analyses will be conducted on the sun-subtracted data. The plots produced will be placed in a subdirectory of
        <input_dir> called "plots". ANYTHING IN THESE DIRECTORIES BEFORE THE SCRIPT IS RUN WILL BE DELETED

    Output: In <input_dir>/plots, there are, as the name suggests, plots. For each of the plots explained below in "Description", there is a .pdf that compiles that plot
        from each pulsar into a single document. In addition, there are some plots that compile the data from all the pulsars, contained in all-pulsar_plots.pdf. Each
        pulsar also has a document that contains all plots relevant to that pulsar. In <input_dir>/dmx, each pulsar has a .dmx file that is the original data with the
        solar wind component subtracted. An additional file is created to store relevant data about all the pulsars called "pulsarStats.txt".

    Description: In order to do the sun-subtraction, the low-frequency components of the dD measurements needs to be removed. This is done by fitting a 3rd-Order Polynomial
        (called f) to the data. In the first set of plots, "poly3fit.pdf", the original dD data (called dD), f, and dD-f are over-plotted. Next, dD-f is folded with a period
        of 1 year to show dD by phase. Next, a model adapted from D. Madison et al. (2018) and J. Hazboun et al. (2021) for the DM contributions due to the solar wind was
        applied (called Z_0). A sinudoid (called S) is then fit to the folded dD-f-Z_0. A plot with folded dD-f, Z_0, S, Z_0+S over-plotted and folded dD-f-Z_0-S with an
        offset for readability is saved in "folded.pdf". Z_0 and S are then unfolded and applied to the complete dD time-series. dD and sun-subtracted dD (dD-Z_0-S) are
        over-plotted in "sun-subtracted_and_original_comparison.pdf". The sun-subtracted dD is output to a .dmx file. Next, a Lomb-Scargle periodogram is created from 
        both dD and sun-subtracted dD. These are over-plotted in "lomb-scargle.pdf". It is normalized by dividing by the mean power. Now, the residuals of dD-f-S-Z_0 will be
        analyzed, as these ideally would be zero if the previous analyses were correct. First, the residuals are plotted by angular separation from the sun for each pulsar in
        "resid_by_ang_sep.pdf". They are then plotted by time in "resid_by_time_series.pdf". Next is the ACF and SF of dD-Z_0-S. These would usually be calculated by shifting the time
        by some time lag, but this only works as expected when the data is sampled at a constant frequency. Thus, an adjustment needs to be made in order to create the ACF
        and SF. Instead, the time-difference between every data point and every other data point (including itself, but no duplicate pairs) is used as the time lag. These time lags
        are binned. Using this technique, the ACF and SF are plotted in "acf.pdf" and "sf.pdf". However, some of the pulsars are very new, and as such do not have many observations
        and therefore time lags. Because of small-number statistics, bins with few points do not have significant data, and are thus removed (in this case, each bin needs at least
        50 points). Additionally, if there aren't enough bins left for the ACF and SF to be meaningful, the plots are not produced (in this case, at least 10 bins). For the structure
        function, it is expected to see a power-law with a slope of 5/3, assuming that the ISM can be treated as a Kolmogorov medium. A powerlaw (a*np.power(x, b)+c) is fit
        to the SF, and a vertical offset is also fit to represent the average noise. This fit is overplotted with the SF. Next, a plot is made to visually show the path of the pulsar
        relative to the sun as a function of time. This plot is saved in "sky_map.pdf". The difference between the narrowband and wideband time series for each pulsar is then
        examined. First, dD_nb, dD_wb, and their difference are over-plotted in "nb-wb_vs_time.pdf". The difference is also plotted as a function of DM, and is saved as "nb-wb_vs_angle.pdf".
        Finally, there are two plots that have data from all pulsars. A version of "resid_by_ang_sep.pdf" is created except with all pulsars over-plotted. This is then transformed
        To a visual map of where the pulsar was measured relative to the sun and what its residual was at that point. These two plots are saved in "all-pulsar_plots.pdf".
        Besides plots, "pulsarStats.txt" is also produced, which acts as a table containing interesting information about each pulsar. The data included are the name of the pulsar,
        whether the specific file was wideband or narrowband, its closest approch to the sun in degrees, the parameters for the fitting functions (3rd-Order Polynomial, Sinusoid,
        and powerlaw for SF), their uncertainties, the weighted RMS deviation of the difference between the wideband and narrowband time series, and the WRMSE of the dD-f-S-Z_0 residuals.
    
    Possible Issues: For f, many of the underlying trends of dD do not follow a 3rd-Order Polynomial, so there might be low-frequency components that are not removed by f subtracting f.
        For the Lomb-Scargle Periodogram, while the power at 1/yr (and harmonics) should decrease, sometimes it does not. This may be due to slight differences in how the two datasets differ, but
        it is probably due to subtracting the solar wind function -- that process not being exact enough.
        Many of the ACFs and SFs look really strange, and additionally many of the fits act strangely. Be cautious when using these results.
        For the sun-subtracted data, some pulsars look much better, and in that case the data produced by this program is preferable to the raw data from the .par file. However, many
        of the pulsars see little changes, and a few may even look (slightly) worse than before. Discretion should be used when using these .dmx files.

    
--solarWindTools.py--
    Description: Contains various function definitions for use in solarWindPlots.py. These include the fitting functions, the solar wind model, and a function to export data as a .dmx file.

--parReader.py--
    Description: A script that converts a .par file into a custom Pulsar Class, which makes usage easier. This was created durint the summer program. It is likely not 100% bug-free,
        optimized, or even useful, but it served its purpose this summer. Perhaps it may be useful to you. The benefit that it has is that in only reads in the data that is available,
        and if the data are missing for certain fields, they are skipped.

--parToDmx.py--
    Description: A script that converts a .par file into a .dmx file. I wasn't sure what to put for average error, so that field has been left blank. Also, it was based off a single .dmx file that
        I referenced, so it may be a slightly different format. This is not actually used in solarWindPlots.py, but I have provided it in hopes that it will be useful.

    Command Line: python3 parToDmx.py <input_name.par> <output_name.dmx>

    Notes: If no output name is given, the input filename will be used.

--solarWindPlots.ipynb--
    Description: A notebook version of solarWindPlots.py that analyzes a single .par file. Used as a way to examine what solarWindPlots.py does.

