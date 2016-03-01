function str=get_strID(cellarr,type)
str='';
switch type
  case {'mech','mechanism','m'}
    mechanisms=cellarr;
    if ischar(mechanisms)
      % convert string '{M1,M2,...}' into cell array {'M1','M2',...}
      mechanisms=regexp(mechanisms,'(\w+)','tokens');
      mechanisms=[mechanisms{:}];
    end
    % use unique to remove duplicates and sort elements
    mechanisms=unique(mechanisms);
    % combine mechanisms into hyphen-separated list 'M1-M2-...'
    for i=1:length(mechanisms)
      str=[str mechanisms{i} '-'];
    end
    str=str(1:end-1);
  case {'vary','v'}
    vary=cellarr;
    % sort vary based on parameter name
    if size(vary,2)==2
      [~,I]=sort(vary(:,1));
    elseif size(vary,2)==3
      [~,I]=sort(vary(:,2));
    end
    vary=vary(I,:);
    % concatenate vary info
    for i=1:size(vary,1)
      if size(vary,2)==2
        str=[str vary{i,1} num2str(vary{i,2}(1)) '-' num2str(vary{i,2}(end)) '_'];
      elseif size(vary,2)==3
        if ~isempty(vary{i,1})
          str=[str vary{i,1} '-' vary{i,2} num2str(vary{i,3}(1)) '-' num2str(vary{i,3}(end)) '_'];
        else
          str=[str vary{i,2} num2str(vary{i,3}(1)) '-' num2str(vary{i,3}(end)) '_'];          
        end
      end
    end    
    str=str(1:end-1);
  case {'mods','modifications'}
    mods=cellarr;
    % sort modifications based on parameter names
    if size(mods,2)==2
      [~,I]=sort(mods(:,1));
    elseif size(mods,2)==3
      [~,I]=sort(mods(:,2));
    end
    mods=mods(I,:);
    % concatenate mod info
    for i=1:size(mods,1)
      if size(mods,2)==2
        str=[str mods{i,1} num2str(mods{i,2}) '_'];
      elseif size(mods,2)==3
        if ~isempty(mods{i,1})
          str=[str mods{i,1} '-' mods{i,2} num2str(mods{i,3}) '_'];
        else
          str=[str mods{i,2} num2str(mods{i,3}) '_'];          
        end
      end
    end
    str=str(1:end-1);
  case {'options'}
    options=cellarr; % key/value options
    for i=1:2:length(options)
      key=options{i}; val=options{i+1};
      if ischar(val)
        str=[str key '-' val '_'];
      elseif isnumeric(val) && length(val)==1
        str=[str key num2str(val) '_'];
      elseif isnumeric(val) && length(val)>1
        str=[str key num2str(val(1)) '-' num2str(val(end)) '_'];
      else
        % don't know what to do...
      end
    end
    str=str(1:end-1);
end


