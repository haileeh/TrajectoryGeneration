# TrajectoryGeneration

This folder contains original optimization code and my re-creation, as well. 

To quickly analyze a dataset that has already been created, run startup, then open Optimization > analyzeOptimized and run.

The main thing to note at this point is that when the Chaser is approaching in the x-direction and the Target is spinning about z (flat spin), the Coriolis acceleration will only be nonzero in y while the centripetal acceleration will only be nonzero in x (and linear acceleration is only nonzero in x). Therefore, the radial acceleration does not equal the tangential acceleration.
