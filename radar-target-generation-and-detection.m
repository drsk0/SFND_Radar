clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
Rmax = 200.0
% Range resolution
Rres = 1.0
% Max Velocity = 70 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8
%% User Defined Range and Velocity of target
R0 = 100.0
v0 = -20.0
% define the target's initial position and velocity. Note : Velocity
% remains contant
 
%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

B = c / (2.0 * Rres)
Tchirp = 5.5 * 2.0 * Rmax/c
alpha = B / Tchirp


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    t_i = t(i);
    R = R0 + v0*t_i; 
    tau = 2*R/c;
    r_t(i) = R;
    td(i) = tau;
    Tx(i) = cos(2 * pi * (fc * t_i + alpha * t_i^2 / 2));
    
    Rx(i)  = cos(2 * pi * (fc * (t_i - tau)) + alpha * (t_i - tau)^2 / 2);
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = cos(4*pi*t_i*(alpha*R + fc*v0) / c);
    
end

%% RANGE MEASUREMENT


%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = Mix + 0.1*randn(size(t));
Mix = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix, Nr, 1);

% Take the absolute value of normalized FFT output
Mix_fft = abs(Mix_fft / Nr);


% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix_fft = Mix_fft(1:Nr/2);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % plot FFT output 
plot(Mix_fft)

 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift(sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 3
Td = 3

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 1
Gd = 1

% offset the threshold by SNR value in dB
offset = 6.0

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2,Nd);
signal_cfar = zeros(Nr/2,Nd);


%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0. 


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR

Nt = (Tr + Gr + 1) * (Td + Gd + 1) - (Gr + 1)*(Gd + 1);
for i=1+Tr+Gr:Nr/2-Tr-Gr
 for j=1+Td+Gd:Nd-Td-Gd
  nl = 0.0;
  for x=i-Tr-Gr:i+Tr+Gr
   for y=j-Td-Gd:j+Td+Gd
    nl = nl + db2pow(RDM(x, y));
   end
  end
  for x=i-Gr:i+Gr
   for y=j-Gd:j+Gd
    nl = nl - db2pow(RDM(x,y));
   end
  end
  nl = nl/Nt;
  nl = pow2db(nl) + offset;
  noise_level(i,j) = nl;
  if RDM(i,j) > nl 
   signal_cfar(i,j) = 1.0;
  end
 end
end





% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 








%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;
