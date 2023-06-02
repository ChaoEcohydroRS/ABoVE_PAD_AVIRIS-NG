
# Mapping vegetation communities in the Arctic-boreal wetlands of the Peace-Athabasca Delta using AVIRIS-NG hyperspectral data

## (Under construction)

Repository named **ABoVE_PAD_AVIRIS-NG** for the paper titled **Mapping vegetation communities in the Arctic-boreal wetlands of the Peace-Athabasca Delta using AVIRIS-NG hyperspectral data** subimitted to RSE.

**Author List**
Chao Wang<sup>1</sup>, Tamlin M. Pavelsky<sup>1</sup>, E. D. Kyzivat<sup>2</sup>, Fenix Garcia-Tigreros<sup>3</sup>, Erika Podest<sup>4</sup>, Fangfang Yao<sup>5</sup>, Xiao Yang<sup>1</sup>, Shuai Zhang<sup>6</sup>, Conghe Song<sup>7</sup>, Theodore Langhorst<sup>1</sup>, Wayana Dolan<sup>1</sup>, Martin R. Kurek<sup>8</sup>, Merritt E. Harlan<sup>9</sup>, Laurence C. Smith<sup>2</sup>, David E. Butman<sup>3</sup>, Robert G. M. Spencer<sup>8</sup>, Colin J. Gleason<sup>9</sup>, Kimberly Wickland<sup>10</sup>, Robert G. Striegl<sup>10</sup>,  Daniel L. Peters<sup>11</sup>

</br>1 Department of Earth, Marine and Environmental Sciences, University of North Carolina, Chapel Hill, NC, USA
</br>2 Department of Earth, Environmental & Planetary Sciences and Institute at Brown for Environment & Society, Brown University, Providence, RI, USA
</br>3 School of Environmental and Forest Sciences, University of Washington, Seattle, WA, USA 
</br>4 Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA
</br>5 Cooperative Institute for Research in Environmental Sciences (CIRES), University of Colorado Boulder, Boulder, CO, USA
</br>6 College of Marine Science, University of South Florida, St. Petersburg, FL, USA 
</br>7 Department of Geography, University of North Carolina, Chapel Hill, NC, USA
</br>8 Earth, Ocean & Atmospheric Science Department, Florida State University, Tallahassee, FL, USA
</br>9 Department of Civil and Environmental Engineering, University of Massachusetts, Amherst, MA, USA
</br>10 Water Resources Mission Area, U.S. Geological Survey, Boulder, CO, USA
</br>11 Department of Geography, University of Victoria, Victoria, BC, Canada



## Abstract
Since climate warming is twice the global average, the northernmost region of the world is suffering from a series of hydrologic changes, such as the permafrost thaw, spring floods, and summer droughts. Studies have shown that unprecedented climate warming has led to shifts in plants’ composition in the boreal ecosystems. It is particularly prominent in the changes of wetland vegetation communities, which affects wetland productivity, species diversity, and biogeochemical cycles. Establishing the relationship between the biophysical characteristics of vegetation communities and the hyperspectral signal provides the possibility of accurately identifying vegetation communities, and thus can complete the landscape-scale analysis of long-term climate-driven vegetation communities responses. However, the spectral signals and characteristics of different wetland vegetations in the boreal ecosystems have not yet been fully explored. Here, based on the manual interpretations of our field drone images and photos in the summer of 2019, we intend to address the applicability of Airborne Visible Infrared Imaging Spectrometer-Next Generation (AVIRIS-NG) hyperspectral data for supervised classifying vegetation communities along the gradients from aquatic-vegetation, emerging-vegetation, to highland-vegetation around the Peace-Athabasca Delta. Our results concluded that the AVIRIS-NG data owing to its high spatial and spectral resolution, yield information on structural details and canopy parameters which has the advanced potential for accurately mapping of vegetation communities of wetland in boreal ecosystems.



# Python Instructions
For this pipeline to work you will need to [sign up](https://earthengine.google.com/signup/) and have a Google Earth Engine configured python installation ready to go. 

<br>Explaining exactly how to do this is beyond the scope of this package but Google provides detailed installation instructions [here](https://developers.google.com/earth-engine/python_install).

## Installation

After installing all necessary package, please check by running the following code:
```python
import os
import sys
import ee
import math
import re
```

### Manually



## Usage
Please email me at chao.wang@unc.edu for any further information.

The GEE repository includes some classification mapping demos:
```javascript
// Store in a variable (see below).

```

## Other Resources
The **Repository** for processing **UAVSAR imagery** can be found in
https://github.com/waynechao128/FlorenceFlood_UAVSAR_Repo
and detailed about this refers to the published paper: Wang, C., Pavelsky, T. M., Yao, F., Yang, X., Zhang, S., Chapman, B., et al. (2022). Flood extent mapping during Hurricane Florence with repeat-pass L-band UAVSAR images. Water Resources Research, 58, e2021WR030606. https://doi.org/10.1029/2021WR030606

## Citation
Wang, C., T. M. Pavelsky, E. D. Kyzivat, F. Garcia-Tigreros, E. Podest, F. Yao, X. Yang, S. Zhang, C. Song, T. Langhorst, W. Dolan, M. R. Kurek, M. E. Harlan, L. C. Smith, D. E. Butman, R. G. M. Spencer, C. J. Gleason, K. P. Wickland, R. G. Striegl, and D. L. Peters. 2023. Quantification of wetland vegetation communities features with airborne AVIRIS-NG, UAVSAR, and UAV LiDAR data in Peace-Athabasca Delta. Remote Sensing of Environment 294:113646. https://doi.org/10.1016/j.rse.2023.113646

## License
The material is made available under the **MIT License**: Copyright 2022, Chao Wang, Tamlin M. Pavelsky, of Global Hydrology Lab - University of North Carolina, Chapel Hill.
All rights reserved.
