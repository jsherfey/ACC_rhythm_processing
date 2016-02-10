function data=DownsampleData(data,downsample_factor)
% post-process downsampling
if downsample_factor>1
  for i=1:length(data)
    for j=1:length(data(i).labels)
      data(i).(data(i).labels{j})=data(i).(data(i).labels{j})(1:downsample_factor:end,:);
      data(i).(data(i).labels{j})=single(data(i).(data(i).labels{j}));
    end
  end
end
