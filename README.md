# Travelling-Salesman-with-Constraint
Travelling Salesman Problem with 30 locations but some locations blocked due to circular storms. Find the best path using optimization.
A cargo truck has to visit 30 locations on an island to deliver cargoes. There is only a single
two-way road to travel from any of these two locations. These roads have no turning points
i.e. they are linear. Suppose that the locations are on a tropical island where regional storms
occur in circular shapes and surround certain regions of the island. If some part of the road is
surrounded by the storm, the cargo truck prefers not to choose that road due to its possible
dangers. Suppose that cargo truck is notified by the exact locations of these storms before it
starts to deliver the cargoes. Moreover, suppose that the cargo truck moves with a constant
speed along the road. There are three types of roads on this island: asphalt, concrete and
gravel. Speed of the truck differs based on the road type:Asphalt Concrete Gravel
100 unit/hour 65 unit/hour 35 unit/hour
you are given the Cartesian coordinates of 30 locations (first sheet of
data.xlsx) and the storm regions defined by the center and radius of the storm (second
sheet of data.xlsx). The type of the road between any of the two locations are also provided
(the third sheet of data.xlsx).
Suppose that the cargo truck drops the cargoes immediately and leaves the location. Hence,
no time delay occurs. The cargo truck’s aim is to visit all the locations from a starting point
(the first location in the first sheet of data.xlsx) and then turn back to its initial location
while minimizing the traveling time.
