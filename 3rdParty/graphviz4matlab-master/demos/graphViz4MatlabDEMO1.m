load smallExample 
nodeColors = {'g','b','r','c'}; % if too few specified, it will cycle through
edgeColors = {'Tom', 'Bill', 'r'
              'Bill' 'all' , 'g'};

graphViz4Matlab('-adjMat',A,'-nodeLabels',T.name_AAL2,'-layout',Treelayout,'-nodeColors',nodeColors,'-edgeColors', edgeColors);
