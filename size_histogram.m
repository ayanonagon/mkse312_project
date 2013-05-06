function [] = size_histogram(modules)
%SIZE_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

bucket = [];
for i=1:length(modules)
    el = modules(i);
    el = el{1};
    bucket = [bucket length(el)];
end

hist(bucket);

