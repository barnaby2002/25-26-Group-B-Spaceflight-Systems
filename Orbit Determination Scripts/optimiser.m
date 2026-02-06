% Constellation sweep: altitude fixed, sweep inclination, total sats, planes, and relative spacing
% MATLAB R2025b + Satellite Communications Toolbox (+ Aerospace Toolbox for walkerDelta)
%
% Goal:
%   Find minimum satellites that achieve continuous >=4-in-view over British Isles
%   with elevation mask = 10 deg, without increasing altitude toward MEO.
%
% What it does:
%   - Builds satelliteScenario once per design point (to avoid stale objects)
%   - Builds Walker-Delta constellation via walkerDelta
%   - Evaluates coverage over a UK grid by counting sats above elevation mask
%   - Records: minInView (global worst-case) and fracGe4Worst (worst gridpoint)
%   - Writes a results table to CSV and prints best candidates.
%
% Tips:
%   - Start coarse, then narrow ranges around good performers.
%   - dt=60s is usually fine for ranking; confirm finalists at dt=10-30s.

clear; clc;

%% Fixed settings
maskDeg   = 10;
alt_km    = 1200;          % keep "low-ish" (edit)
dt        = 60;            % [s] coarse sweep; tighten later
daysSim   = 1;             % simulate 1 day
startTime = datetime(2035,1,1,0,0,0,"TimeZone","UTC");
stopTime  = startTime + days(daysSim);

Re = 6378.137e3;
a  = Re + alt_km*1e3;
e  = 0;

%% UK grid (coarse). Refine later.
latVec = 49:1:61;
lonVec = -8:1:2;
[Lon, Lat] = meshgrid(lonVec, latVec);
nPts = numel(Lat);

%% Sweep ranges (edit these)
incList = 60:2:78;                 % deg
tList   = [60 72 84 88 92 96];     % total satellites
pList   = 3:12;                    % planes (only those dividing t will be used)
fList   = 0:5;                     % relative spacing candidates (filtered to < p later)

%% Pass/fail requirement
reqMinInView = 4;                  % continuous >=4 (global min)
reqFracWorst = 1.0;                % worst point fraction >=4 (continuous)

%% Preallocate results (worst-case size)
maxRows = numel(incList)*numel(tList)*numel(pList)*numel(fList);
R = table('Size',[0 8], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'alt_km','inc_deg','t','p','f','minInView','fracGe4Worst','fracGe4Median'});

row = 0;

%% Sweep
for inc = incList
    for t = tList
        for p = pList
            if mod(t,p) ~= 0
                continue;
            end

            % f is usually meaningful in [0, p-1]
            for f = fList
                if f >= p
                    continue;
                end

                % Fresh scenario each design point (avoids reusing old sats/gs)
                sc = satelliteScenario(startTime, stopTime, dt);

                % Constellation
                sats = walkerDelta(sc, a, inc, t, p, f);
                Ns = numel(sats);

                % Ground stations: create per-point in this scenario
                gs = cell(nPts,1);
                for q = 1:nPts
                    gs{q} = groundStation(sc, ...
                        "Name", sprintf("GP_%03d", q), ...
                        "MinElevationAngle", maskDeg, ...
                        "Latitude", Lat(q), ...
                        "Longitude", Lon(q));
                end

                % Sample times (scenario cadence)
                tvec = (sc.StartTime:seconds(dt):sc.StopTime).';
                Nt = numel(tvec);

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

                minInView_global   = min(minInView_pts);
                fracGe4Worst       = min(fracGe4_pts);
                fracGe4Median      = median(fracGe4_pts);

                row = row + 1;
                R(row,:) = {alt_km, inc, t, p, f, minInView_global, fracGe4Worst, fracGe4Median};

                fprintf("alt=%4g km inc=%5.1f t=%3d p=%2d f=%d  -> minInView=%d  worstFrac>=4=%.3f\n", ...
                    alt_km, inc, t, p, f, minInView_global, fracGe4Worst);
            end
        end
    end
end

%% Filter passing designs
pass = R.minInView >= reqMinInView & abs(R.fracGe4Worst - reqFracWorst) < 1e-12;
Rpass = R(pass,:);

% Sort: minimum satellites first, then fewer planes, then higher worst fraction (tie-break)
if ~isempty(Rpass)
    Rpass = sortrows(Rpass, {'t','p','inc_deg','f'});
    disp("=== Passing designs (continuous >=4 everywhere on grid) ===");
    disp(Rpass(1:min(20,height(Rpass)),:));
else
    disp("No passing designs found in this sweep range. Expand tList and/or alt_km, or reduce dt/grid coarseness.");
end

%% Save results
outCsv = sprintf("sweep_results_alt%dkm_mask%d_dt%ds.csv", alt_km, maskDeg, dt);
writetable(R, outCsv);
fprintf("Wrote %d rows to %s\n", height(R), outCsv);

%% Optional: quick "best near-pass" list
% Sort by minInView desc, worstFrac desc, then satellites asc
Rnear = sortrows(R, {'minInView','fracGe4Worst','t','p'}, {'descend','descend','ascend','ascend'});
disp("=== Best near-pass candidates (top 20) ===");
disp(Rnear(1:min(20,height(Rnear)),:));
