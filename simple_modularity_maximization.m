function [groups] = simple_modularity_maximization(adj)
%SIMPLE_MODULARITY_MAXIMIZATION Summary of this function goes here
%   Detailed explanation goes here

groups={};
currentOverallBest = -1;
previousOverallBest = -1;
bestMove = -1;
bestdQ = -Inf;
n=size(adj,1); % number of vertices
Qs = [];
states = [];

for i=1:n
  grp = mod(i, 2) + 1;
  if length(groups) > 1
      groups{grp}=[groups{grp} i];
  else
     groups{grp}=[i]; 
  end
end

currentGroups = groups;

currentOverallBest = modularity_metric(groups, adj);

moved = [];
while currentOverallBest ~= previousOverallBest
    if length(moved) >= size(adj, 1)
        moved = [];
        [~, idx] = max(Q);
        currentGroups = states(idx);
    end
    % find best to move
    Q = modularity_metric(currentGroups, adj);
    for i=1:n
        if isempty(find(moved == i, 1))
            disp(i);
            newQ = modularity_metric(move_item(currentGroups, i), adj);
            newdQ = newQ - Q;
            if newQ - Q > bestdQ
                bestMove = i;
                bestdQ = newdQ;
            end
        end
    end
    currentGroups = move_item(currentGroups, bestMove);
    states = [states currentGroups];
    currentQ = modularity_metric(groups, adj);
    Qs = [Qs currentQ];
    previousBest = currentOverallBest;
    currentOverallBest = currentQ;
    disp(currentOverallBest - previousBest);
end
end

function [groups] = move_item(groups, i)
if which_group(groups, i) == 1
    g1 = groups{1};
    groups{1} = g1(g1~=i);
    groups{2} = [groups{2} i];
else
    g2 = groups{2};
    groups{2} = g2(g2~=i);
    groups{1} = [groups{1} i];
end

end

function [groupNum] = which_group(groups, i)
if ~isempty(find(groups{1} == i, 1))
    groupNum = 1;
else
    groupNum = 2;
end
end