
%reading data from the files
location_set=xlsread('data.xlsx','Locations','B2:C31'); %read location coordinates
storm_set = xlsread('data.xlsx','Storms','A3:C22'); %read storm locations
[road_paths,Road_type,raw_road_data] = xlsread('data.xlsx', 'Road Material','A2:C871'); %read road data
stopsLon = location_set(:,1); %all x coordinates of the location
stopsLat=location_set(:,2); %all y coordinates of the location
x_circle=storm_set(:,1); %all x coordinates of the storms
y_circle=storm_set(:,2); %all y coordinates of the storms
radius_circle=storm_set(:,3); %all radius of the storms
centress = [x_circle y_circle]; %the coordinates of the storms

nStops=30; %number of locations

% plotting Circular Storms and all Locations
hold on
axis equal
figure(1)
scatter(stopsLon,stopsLat, 'filled') %plot all locations
scatter(167.2313144,-112.331426,'g','filled') %distinguish starting location in different color
viscircles(centress,radius_circle) %plot all circular storms
title('Delivery Points and Circular Storms')
%legend('Locations','Starting Location')
hold off

Possible_Trips = nchoosek(1:nStops,2); %generate all trips / all pairs of stops
r=size(Possible_Trips,1);

%identify the paths which have storms
i=size(location_set,1); %size of the location matrix
u=size(storm_set,1); %size of the storm matrix

A_path= [];
B_path=[];

%using intersection to find intersection points between the paths and the
%circles
for x= 1:i
    for y= 1:i
        coord1 = location_set(x,:);
        coord2 = location_set(y,:);
        xdata1 = [coord1(1) coord2(1)];
        ydata1 = [coord1(2) coord2(2)];
        t=linspace(0,2*pi,30);
        for j = 1:u
            x2=radius_circle(j)*cos(t) + x_circle(j);
            y2=radius_circle(j)*sin(t)+y_circle(j);
            [xc,yc] = polyxpoly(xdata1,ydata1,x2,y2);
            if [xc,yc]~=0
                A_path= [A_path;{coord1}];
                B_path=[B_path;{coord2}];
                break
            end
        end
    end
end

s=size(A_path,1);

A_path_node=[];
%identifying paths which have storms by node number
for tyo= 1:s
    for tyu = 1:i
        if A_path{tyo}==location_set(tyu,:)
            A_path_node = [A_path_node;tyu];
        end
    end
end

B_path_node=[];
for gu= 1:s
    for cb = 1:i
        if B_path{gu}==location_set(cb,:)
            B_path_node = [B_path_node;cb];
        end
    end
end

combined_nodes_paths= cat(2,A_path_node,B_path_node);

refined_combined_nodes_paths= [];

for wt = 1:s
    for ht = 1:r
        if Possible_Trips(ht,:)==combined_nodes_paths(wt,:)
            refined_combined_nodes_paths = [refined_combined_nodes_paths;combined_nodes_paths(wt,:)];
        end
    end
end

%find Nodes which do not have Storms
Remaining_Nodes=nchoosek(1:nStops,2);
Remaining_Nodes(ismember(Remaining_Nodes,refined_combined_nodes_paths,'rows'),:)=[];

%Plot all paths which have circular Storms
rty=location_set(A_path_node,:);
truy=location_set(B_path_node,:);
for yyh = 1:s
    hold on
    axis equal
    plot([rty(yyh) truy(yyh)],[rty(yyh,2) truy(yyh,2)])
    figure(2)
    title('All paths which are blocked due to Storm')
end
scatter(stopsLon,stopsLat, 'filled')
scatter(167.2313144,-112.331426,'g','filled')
viscircles(centress,radius_circle)
hold off

%calculate the trip distances
dist = hypot(stopsLat(Remaining_Nodes(:,1)) - stopsLat(Remaining_Nodes(:,2)), ...
    stopsLon(Remaining_Nodes(:,1)) - stopsLon(Remaining_Nodes(:,2)));
lendist = length(dist);

k=size(Remaining_Nodes,1);
b=size(raw_road_data,1);

roadtype_allowed = strings;
raw_road_data =(string(raw_road_data));

for gg=1:k
    for bb=1:b
        if string(Remaining_Nodes(gg,:))==(raw_road_data(bb,1:2))
            roadtype_allowed=cat(1,roadtype_allowed,raw_road_data(bb,3));
        end
    end
end
roadtype_allowed(1) = [];

speeds_of_allowed_nodes = [];

for u = 1:size(roadtype_allowed,1)
    if roadtype_allowed(u)=='Asphalt'
        speeds_of_allowed_nodes = [speeds_of_allowed_nodes ; 100];
    elseif roadtype_allowed(u)=='Concrete'
        speeds_of_allowed_nodes = [speeds_of_allowed_nodes ; 65];
    else
        speeds_of_allowed_nodes = [speeds_of_allowed_nodes ; 35];
    end
end

%calculate time for the allowed paths
time_allowed_paths = dist./speeds_of_allowed_nodes;

%create a graph where stops are nodes and trips are edges
G = graph(Remaining_Nodes(:,1),Remaining_Nodes(:,2));
figure(3)
axis equal
hGraph = plot(G,'XData',stopsLon,'YData',stopsLat,'LineStyle','none','NodeLabel',{})
title('Optimal Path after all Constraints')

% create binary optimization variable representing the potential trips.
tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);

%objective function is
tsp.Objective = time_allowed_paths'*trips;

%add linear constraints that each stop has two assosiated trips i.e. there must be a trip to and from each stop
%identify all trips starting or ending at a stop by finding all edges connecting to that stop
%For each stop, create the constraint that the sum of trips for that stop equals two.
constr2trips = optimconstr(nStops,1);
for stop = 1:nStops
    whichIdxs = outedges(G,stop); % Identify trips associated with the stop
    constr2trips(stop) = sum(trips(whichIdxs)) == 2;
end
tsp.Constraints.constr2trips = constr2trips;

%solve problem
opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts);

%Create a new graph with the solution trips as edges
%round the solution in case some values are not exactly integers, and convert the resulting values to logical
tspsol.trips = logical(round(tspsol.trips));
Gsol = graph(Remaining_Nodes(tspsol.trips,1),Remaining_Nodes(tspsol.trips,2));
%Overlay the new graph on the existing plot and highlight its edges.
hold on
highlight(hGraph,Gsol,'LineStyle','-');

%Subtours Constraints using Iterative approach
% Detect the subtours in the current solution,
%then add inequality constraints to prevent those particular subtours from happening.
tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % Number of subtours
fprintf('# of subtours: %d\n',numtours);
%5 subtours are found

%Include the linear inequality constraints to eliminate subtours,
%and repeatedly call the solver, until just one subtour remains

% Index of added constraints for subtours
kbnm = 1;
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    for ii = 1:numtours
        inSubTour = (tourIdxs == ii); % Edges in current subtour
        a = all(inSubTour(Remaining_Nodes),2); % Complete graph indices with both ends in subtour
        constrname = "subtourconstr" + num2str(kbnm);
        tsp.Constraints.(constrname) = sum(trips(a)) <= (nnz(inSubTour) - 1);
        kbnm = kbnm + 1;
    end
    
    % Try to optimize again
    [tspsol,fval,exitflag,output] = solve(tsp,'options',opts);
    tspsol.trips = logical(round(tspsol.trips));
    Gsol = graph(Remaining_Nodes(tspsol.trips,1),Remaining_Nodes(tspsol.trips,2));
    
    % Plot new solution
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-');
    axis equal
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % Number of subtours
    fprintf('# of subtours: %d\n',numtours)
end

disp(output.absolutegap)











 
