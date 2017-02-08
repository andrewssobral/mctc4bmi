function perform_frame_selection(params)
  %% Load images
  clc;
  sequence_name = params.sequences_name(params.current_sequence).name;
  disp(['Loading images from sequence: ' sequence_name]);
  V = load_images(params); % show_4dvideo(V);
  disp(size(V));
  
  %%% Frame difference
  disp('Performing frame difference');
  O = zeros(size(V,1),size(V,2),size(V,4));
  vDiffSum = zeros(size(V,4),1);
  myfilter = fspecial('gaussian',[3 3], 0.5);
  I0 = rgb2gray(V(:,:,:,1));
  switch(sequence_name)
    case 'Toscana'
      beta = 1e-2;
    otherwise
      beta = 1e-3;
  end
  for i = 2 : size(V,4)
    if(mod(i,10) == 0) fprintf('%d ',i); end
    
    I = rgb2gray(V(:,:,:,i));
    I = imfilter(I, myfilter, 'replicate');
    %I = medfilt2(I,[5,5]);
    %I = imgaussfilt(I,2);
    
    switch(sequence_name)
      case {'Candela_m1.10','HallAndMonitor'}
        if(i == 2) disp('Using first image as reference...'); end
        I_old = I0;
      otherwise
        if(i == 2) disp('Using the previous frame as reference...'); end
        I_old = rgb2gray(V(:,:,:,i-1));
    end
    I_old = imfilter(I_old, myfilter, 'replicate');

    I_mag = sqrt(power(I-I_old,2));
    I_res = I_mag;

    I_thr = 0.5*I_res.^2 > beta;
    O(:,:,i) = I_thr;
    vDiffSum(i) = sum(I_thr(:));
  end
  disp(i);
  vDiffSum = vDiffSum./max(vDiffSum(:));
  % show_3dvideo(O);
  clear i I I0 I_old I_mag I_res I_thr myfilter beta;
  
  %%% Average
%   switch(sequence_name)
%     case 'HallAndMonitor'
%       avgDiffSum = medfilt1(vDiffSum,5,'zeropad');
%     otherwise
       avgDiffSum = medfilt1(vDiffSum,3,'zeropad');
%   end
  if(params.debug)
    clf,plot(vDiffSum),hold on,plot(avgDiffSum),hold off;
    pause(1);
  end
  
  %% Gradient
  disp('Calculating difference gradient');

  %S = vDiffSum;
  S = avgDiffSum;
  
  Y = gradient(S);
  %Y = diff(S);
  Y = (Y-min(Y))/(max(Y)-min(Y)); % normalize
  %Y = medfilt1(Y,3,'truncate');
  %Y = medfilt1(Y,3,'zeropad');

  switch(sequence_name)
    case {'HallAndMonitor','HighwayII'}
      Tt = 0.075;
    case {'CAVIAR1','Candela_m1.10','HighwayI','IBMtest2','PeopleAndFoliage'}
      Tt = 0.1;
    case {'Foliage'}
      Tt = 0.2;
    case {'HumanBody2'}
      Tt = 0.05;
    case {'Board','CAVIAR2','CaVignal','Snellen','Toscana'}
      Tt = 0.125;
  end

  T1 = (Y > (mean(Y) + Tt));
  T2 = (Y < (mean(Y) - Tt));
  T = double(or(T1,T2));

  if(params.debug)
    clf;
    plot(S);
    hold on;
    plot(Y);
    plot(T.*.05,'s');
    hold off;
    title('Frame Selection');
    xlabel('Frames'),ylabel({'Difference between';'consecutive frames'});
    legend('vector d normalized ','derivative of vector d',...
      'selected frames',...
      'Location','southoutside',... %northoutside southoutside
      'Orientation','horizontal');
    pause(1);
  end
  disp(['Number of selected frames: ' num2str(sum(T))]);
  
  %% Frame selection
  disp('Performing frame selection');
  switch(sequence_name)
    case {'Toscana'}
      V2 = V;
      O2 = O;
    otherwise
      V2 = V(:,:,:,logical(T));
      O2 = O(:,:,logical(T));
  end
  % size(Vnew), size(Onew)
  
  if(params.debug)
    clf,show_4dvideo(V2);
    %clf,show_3dvideo(O2);slice3(O2);
  end
  
  %% Save
  clear T Tt S T1 T2 Y avgDiffSum;
  disp(['Saving at: ' fullfile(params.sequences_mat,[sequence_name '.mat'])]);
  save(fullfile(params.sequences_mat,[sequence_name '.mat']));
end
