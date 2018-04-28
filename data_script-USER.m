%%
names = dir('*/*.wav');
siz=length(names);
wind = cell(1,siz);
for i=1:length(names)
    wind(i) = {audioread(names(i).name)};
end
%{ 
Rondane 1 and 2 clean wind noises- from rondane which has speech.
 - 00:12 and 00:22 respectively. 
bare poles - wind with wave noise - 1:47
trees - wind with leaves sounds - 00:39
mountain top - wind with a buzz - 00:44
a cut from sous le soleil - wind with ocean sounds and some speech - 00:34
%}
%%
k = 6;
fs=44.1e3;
siz1=length(wind{k});
time=(0:siz1-1)./fs;
plot(time,wind{k}(:,1))
%soundsc(wind{1},fs)
%clear sound