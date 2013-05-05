load 'facebook_adj.mat';

seenA = [];
for a=1:length(A)
    disp(a);
    if ~isempty(find(seenA == a))
        continue;
    end
    seen = [];
    q = [a];
while ~isempty(q) && length(seen) < 2500
    i = q(1);
    q = q(2:end);
    for j = find(A(i, :))
        if isempty(find(seen == j, 1))
            q = [q j];
        end
    end
    seen = [seen i];
end

if length(seen) >= 2500
    disp(seen);
end

seenA = [seenA seen];
disp(seenA);

end