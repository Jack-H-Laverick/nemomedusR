# nemomedusR 1.1.0.0 
<span style="color:grey;">28/11/2022</span>

* Added scheme_column.
* Added "collect" analysis to NEMO_MEDUSA() for subsetting tasks. 

# nemomedusR 1.0.0.1 
<span style="color:grey;">10/02/2021</span>

* Updated tests to use new column names.

# nemomedusR 1.0.0.0 
<span style="color:grey;">02/02/2021</span>

<strong> Major overhaul with breaking changes <strong/>

* The back end has been rewritten using Rcpp for a speed boost and a streamlined codeset to maintain. See the slabR package.
* The summary routines are now accessed through the function NEMO_MEDUSA(), not whole_month().
* slabR has superseded the space and crop arguments with the use of summary "schemes".
* A summary scheme helper has been created for StrathE2E model summaries scheme_strathE2E().
* A summary scheme helper has been created for linear interpolation between depth layers scheme_interp_slice().
* type_in_month() has been replaced by temporal_operation() to more accurately reflect what it does and open up a general implementation for summarising (or not) time steps.
* Extractor functions are now exported functions so users can see what's going on under the hood of NEMO_MEDUSA().
* The documentation for extractors has been streamlined by creating a shared help page per analysis type.
* A number of functions which were only tangentially related to NEMO-MEDUSA model output have been moved back to MiMeMo.tools.

# nemomedusR 0.1.0.9000 
<span style="color:grey;">26/12/2020</span>

* Initialised a test suite.

# nemomedusR 0.0.0.9000 
<span style="color:grey;">09/10/2020</span>

* Excised NEMO-MEDUSA functions from *MiMeMo.tools*.
