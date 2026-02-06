%% Fixed settings
maskDeg   = 10;
dt        = 60;                 % [s] coarse sweep
daysSim   = 100;
startTime = datetime(2035,1,1,0,0,0,"TimeZone","UTC");
stopTime  = startTime + days(daysSim);
endTime = datetime(2040,1,1,0,0,0,"TimeZone","UTC");

Re = 6378.137e3;                % [m]

% UK grid
latVec = 49:1:61;
lonVec = -8:1:2;
[Lon, Lat] = meshgrid(lonVec, latVec);
nPts = numel(Lat);

% requirement
reqMinInView = 4;

%% Results table
R = table('Size',[0 8], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'alt_km','inc_deg','t','p','f','minInView','fracGe4Worst','fracGe4Median'});
row = 0;

gs = cell(nPts,1);
for q = 1:nPts
    gs{q} = groundStation(sc, ...
        "Name", sprintf("GP_%03d", q), ...
        "MinElevationAngle", maskDeg, ...
        "Latitude", Lat(q), ...
        "Longitude", Lon(q));
end

% Compute metrics per grid point
minInView_pts = zeros(nPts,1);
fracGe4_pts   = zeros(nPts,1);

for q = 1:nPts
    vis = false(Nt, Ns);
    for s = 1:Ns
        [~, el] = aer(gs{q}, sats(s));
        vis(:,s) = (el >= maskDeg);
    end

    inView = sum(vis, 2);
    minInView_pts(q) = min(inView);
    fracGe4_pts(q)   = mean(inView >= reqMinInView);
end

minInView_global = min(minInView_pts);
fracGe4Worst     = min(fracGe4_pts);
fracGe4Median    = median(fracGe4_pts);

row = row + 1;
R(row,:) = {alt_km, inc_deg, t, p, f, minInView_global, fracGe4Worst, fracGe4Median};

fprintf("alt=%4.0f km inc=%5.1f t=%3d p=%2d f=%2d -> minInView=%d worstFrac>=4=%.3f\n", ...
    alt_km, inc_deg, t, p, f, minInView_global, fracGe4Worst);

%% Passing designs (continuous at sample times)
pass = (R.minInView >= reqMinInView) & (abs(R.fracGe4Worst - 1.0) < 1e-12);
Rpass = R(pass,:);

if ~isempty(Rpass)
    Rpass = sortrows(Rpass, {'t','p','alt_km','inc_deg','f'});
    disp("=== Passing designs (continuous >=4 everywhere on grid at sampled times) ===");
    disp(Rpass(1:min(20,height(Rpass)),:));
else
    disp("No passing designs found. Increase t / altitude, broaden planes, or reduce mask.");
end
