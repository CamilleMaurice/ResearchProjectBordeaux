% This works on 64bit linux
mex -largeArrayDims 'maxflowmex.cpp' 'maxflow-v3.0/graph.cpp' 'maxflow-v3.0/maxflow.cpp'
%mex -largeArrayDims mex_maxflow.cpp maxflow.cpp graph.cpp 
