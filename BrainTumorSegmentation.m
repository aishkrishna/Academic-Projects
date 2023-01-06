
clear all;
close all;

load("braintumor2ddata.mat");
path = pwd;
myfolder = 'output_flairall';
folder = mkdir([path,filesep,myfolder]);
path  = [path,filesep,myfolder];


for N= 1:3


%disp(flair_all)
    %N= 6 % Patient number 6
    image = squeeze(flair_all(:,:,N))
    image = imrotate(image, 270)
    %I = imshow(image,[])
    
    figure(N);
    set(gcf, 'Color', [0.5,0.5,0.5], 'position', [434,380,1227,597])
    subplot(2,4,1)
    imshow(imresize(image, [256,256]),[]);
    title('Original image','FontSize',15);

    im = im2double(image)
    whos im
    imwrite(im, '/Users/aishkrishnamohan/Downloads/Project3Data_work/patient6.jpg')
    
    brain = imresize(imread('patient6.jpg'), [256,256])
    
    %Filtering
    K = imdiffusefilt(brain)
    K= uint8(K);
    
    K=imresize(K,[256,256]);
    if size(K,3)>1
        K=rgb2gray(K);
    end
    subplot(2,4,2)
    imshow(K, []);
    title('Filtered image','FontSize',15);
    
    %Thresholding segmentation
    out=imresize(K,[256,256]);
    t0=mean(im(:));
    th=t0+((max(K(:))+min(K(:)))./2);
    for i=1:1:size(K,1)
        for j=1:1:size(K,2)
            if K(i,j)>th
                out(i,j)=1;
            else
                out(i,j)=0;
            end
        end
    end
    
    label=bwlabel(out);
    stats=regionprops(logical(out),'Solidity','Area','BoundingBox');
    density=[stats.Solidity];
    area=[stats.Area];
    high_dense_area=density>0.2;
    max_area=max(area(high_dense_area));
    tumor_label_th=find(area==max_area);
    tumor_th=ismember(label,tumor_label_th);
    
    if max_area>50
       subplot(2,4,3)
       imshow(tumor_th,[])
       title('Detected Tumor from Thresholding','FontSize',15);
    else
        h = msgbox('Tumor is not present!','status');
        disp('no tumor');
    end
    
   
    % Watershed segmentation
    
    D = -bwdist(~out);
    L = watershed(D);
    L(~out) =0;
    rgb1 = label2rgb(L, 'parula')
    
    
    %figure;
    subplot(2,4,5); imshow(K, []);
    title('Filtered Image', 'FontSize', 15)
    P = imfuse(K, rgb1, 'blend', 'Scaling', 'joint')
    %subplot(1,4,2);imshow(D,[]);
    subplot(2,4,6); imshow(P);
    title('Watershed Segmentation', 'FontSize', 15);
    

    SE =  [0 1 1; 1 1 1; 0 1 0];
    J = imerode(L, SE)
    %figure;
    %subplot(2,3,6);
    %imshow(J);

    %Morphological Operations
    
    label=bwlabel(J);
    stats=regionprops(logical(J),'Solidity','Area','BoundingBox');
    density=[stats.Solidity];
    area=[stats.Area];
    high_dense_area=density>0.3;
    max_area=max(area(high_dense_area));
    tumor_label=find(area==max_area);
    tumor=ismember(label,tumor_label);
    
    if max_area>150
       %subplot(1,3,3)
       rgb = label2rgb(tumor,'parula')
       subplot(2,4,7);
       imshow(rgb)
       title('Tumor alone after watershed','FontSize',15);
    else
        h = msgbox('Tumor is not present!','status');
        disp('no tumor');
        
    end
    
    %% Getting Tumor Outline - image filling, eroding, subtracting USING MORPHOLOGICAL OPERATIONS
    % erosion the walls by a few pixels
    
    dilationAmount = 3;
    rad = floor(dilationAmount);
    [r,c] = size(tumor_th);
    filledImage = imfill(tumor_th, 'holes');
    
    for i=1:r
       for j=1:c
           x1=i-rad;
           x2=i+rad;
           y1=j-rad;
           y2=j+rad;
           if x1<1
               x1=1;
           end
           if x2>r
               x2=r;
           end
           if y1<1
               y1=1;
           end
           if y2>c
               y2=c;
           end
           erodedImage(i,j) = min(min(filledImage(x1:x2,y1:y2)));
       end
    end
    %figure
    %imshow(erodedImage);
    %title('eroded image','FontSize',20);
    
    %% subtracting eroded image from original BW image
    
    tumorOutline=tumor_th;
    tumorOutline(erodedImage)=0;
    
    %figure;  
    %imshow(tumorOutline);
    %title('Tumor Outline','FontSize',20);
    
    %% Inserting the outline in filtered image in red color
    
    rgb = K(:,:,[1 1 1]);
    red = rgb(:,:,1);
    red(tumorOutline)=255;
    green = rgb(:,:,2);
    green(tumorOutline)=0;
    blue = rgb(:,:,3);
    blue(tumorOutline)=0;
    
    tumorOutlineInserted(:,:,1) = red; 
    tumorOutlineInserted(:,:,2) = green; 
    tumorOutlineInserted(:,:,3) = blue; 
    
    C = imfuse(P,tumorOutlineInserted,'blend', 'Scaling', 'joint')
    %figure;

    subplot(2,4,4);
    imshow(C);
    
    title('Detected Tumor from Watershed','FontSize',15);

    
    image_t = squeeze(truth_all(:,:,N))
    image_t = imrotate(image_t, 270)
    %I = imshow(image_t,[])
    %im_t = im2double(image_t)
    %whos im_t
    %imwrite(im_t, '/Users/aishkrishnamohan/Downloads/Project3Data_work/patient1.jpg')
    
    %brain_t = imresize(imread('patient1.jpg'), [256,256])
    
    subplot(2,4,8);
    imshow(imresize(image_t, [256,256]),[]);
    title('Truth image','FontSize',15);

    temp=[path,filesep,'fig',num2str(N),'.png'];
    saveas(gca,temp);
end


