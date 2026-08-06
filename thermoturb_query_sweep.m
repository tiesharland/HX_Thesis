function sweep = thermoturb_query_sweep(Re_min, Re_max, Pr, Tb, Tw, flow_model)
% Queries thermoturb.com with generateCoefficients=on to get a full Re
% sweep from Re_min to Re_max at fixed (Pr, Tb, Tw) for the specified
% flow model ('laminar' or 'turbulent').
%
% Returns a table with columns:
%   Re_bulk, Nu, Cf, Stanton, Retau, Retau_cp

    import matlab.net.http.*
    import matlab.net.http.field.*
    import matlab.net.URI

    base_url = 'https://www.thermoturb.com';

    %% ---- Build POST body ---- %%

    body = strjoin([
        "fluidModel=variable"
        "flowModel="         + flow_model
        "flowType=pipe"
        "heatForcingType=UIH"
        "thermalBoundaryCondition=uniformHeating"
        "targetType=1"
        "targetValue="       + num2str(round((Re_min+Re_max)/2), '%.0f')
        "generateCoefficients=on"
        "minReynoldsNumber=" + num2str(Re_min, '%.0f')
        "maxReynoldsNumber=" + num2str(Re_max, '%.0f')
        "PrandtlNumber="     + num2str(Pr,     '%.4f')
        "bulkTemperature="   + num2str(Tb,     '%.2f')
        "wallTemperature="   + num2str(Tw,     '%.2f')
        "fluid=Air"
        "viscosity=Sutherland"
        "dynamicInput=110.4"
    ], '&');

    %% ---- POST /generate_plot ---- %%

    post_headers = [
        HeaderField('Content-Type', 'application/x-www-form-urlencoded')
        HeaderField('Referer', [base_url '/input'])
        HeaderField('Origin',  base_url)
    ].';

    post_req  = RequestMessage('POST', post_headers, char(body));
    post_resp = send(post_req, URI([base_url '/generate_plot']));

    if post_resp.StatusCode ~= 200
        error('thermoturb_query_sweep:postFailed', ...
              'POST /generate_plot failed with status %d', ...
              double(post_resp.StatusCode));
    end

    set_cookie_field = post_resp.getFields('Set-Cookie');
    if isempty(set_cookie_field)
        error('thermoturb_query_sweep:noCookie', ...
              'No Set-Cookie header found in POST response.');
    end
    cookie_value = regexp(set_cookie_field.Value, '^[^;]+', 'match', 'once');

    %% ---- GET /download_data ---- %%

    cb       = num2str(round(posixtime(datetime('now')) * 1000));
    get_req  = RequestMessage('GET', HeaderField('Cookie', cookie_value));
    get_resp = send(get_req, URI([base_url '/download_data?cb=' cb]));

    if get_resp.StatusCode ~= 200
        error('thermoturb_query_sweep:getFailed', ...
              'GET /download_data failed with status %d', ...
              double(get_resp.StatusCode));
    end

    %% ---- Unzip response ---- %%

    zip_bytes = get_resp.Body.Data;
    tmp_zip   = [tempname() '.zip'];
    fid = fopen(tmp_zip, 'wb');
    fwrite(fid, zip_bytes);
    fclose(fid);

    tmp_dir = tempname();
    mkdir(tmp_dir);
    unzip(tmp_zip, tmp_dir);

    %% ---- Parse coefficients.dat ---- %%
    % File structure: alternating header line / data line pairs
    % Columns: Retau, Rebulk, Cf, Stanton, Nusselt, Retau_cp

    dat_path = fullfile(tmp_dir, 'coefficients.dat');
    raw      = readlines(dat_path);
    raw      = raw(strlength(strtrim(raw)) > 0);

    is_header  = contains(raw, 'Retau') & contains(raw, 'Rebulk');
    data_lines = raw(~is_header);

    n        = numel(data_lines);
    Retau    = zeros(n, 1);
    Re_bulk  = zeros(n, 1);
    Cf       = zeros(n, 1);
    Stanton  = zeros(n, 1);
    Nu       = zeros(n, 1);
    Retau_cp = zeros(n, 1);

    for i = 1:n
        vals = sscanf(data_lines(i), '%f');
        if numel(vals) >= 6
            Retau(i)    = vals(1);
            Re_bulk(i)  = vals(2);
            Cf(i)       = vals(3);
            Stanton(i)  = vals(4);
            Nu(i)       = vals(5);
            Retau_cp(i) = vals(6);
        end
    end

    %% ---- Remove failed parses and sort by Re ---- %%

    valid    = Re_bulk > 0;
    Re_bulk  = Re_bulk(valid);
    Cf       = Cf(valid);
    Stanton  = Stanton(valid);
    Nu       = Nu(valid);
    Retau    = Retau(valid);
    Retau_cp = Retau_cp(valid);

    [Re_bulk, idx] = sort(Re_bulk);
    Cf       = Cf(idx);
    Stanton  = Stanton(idx);
    Nu       = Nu(idx);
    Retau    = Retau(idx);
    Retau_cp = Retau_cp(idx);

    if isempty(Re_bulk)
        error('thermoturb_query_sweep:noData', ...
              'No valid data rows found in coefficients.dat — check server response.');
    end

    fprintf('    Parsed %d Re points  (Re=[%.0f - %.0f]  Nu=[%.3f - %.3f])\n', ...
            numel(Re_bulk), min(Re_bulk), max(Re_bulk), min(Nu), max(Nu));

    %% ---- Pack into table ---- %%

    sweep = table(Re_bulk, Nu, Cf, Stanton, Retau, Retau_cp);

    %% ---- Cleanup temp files ---- %%

    delete(tmp_zip);
    rmdir(tmp_dir, 's');

end