# qml-pub-bias
Estimates t-hurdles, bias adjustments, and local FDRs on the Chen-Zimmermann Dataset for the paper "Do t-stat hurdles need to be raised?"

This version uses QML, replacing several previous versions that use SMM.  QML is much cleaner in this setting due to the huge amount of estimation noise that comes with truncated simulations (unless you do the truncated simulations very carefully).

Having said that the results using QML and SMM are not very different.  The most recent SMM version is here https://github.com/chenandrewy/t-hurdles
