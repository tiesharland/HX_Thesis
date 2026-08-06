function [Re_b, Nu_out, Cf_out] = thermoturb_query(Re_bulk, Pr, Tb, Tw)

    import matlab.net.http.*
    import matlab.net.http.field.*
    import matlab.net.URI

    base_url = 'https://www.thermoturb.com';

    %% ---- Build form-urlencoded body ---- %%

    body = strjoin([
        "fluidModel=variable"
        "flowModel=laminar"
        "flowType=pipe"
        "heatForcingType=UIH"
        "thermalBoundaryCondition=uniformHeating"
        "targetType=1"
        "targetValue=" + num2str(Re_bulk,'%.6f')
        "minReynoldsNumber=1"
        "maxReynoldsNumber=5000"
        "PrandtlNumber=" + num2str(Pr,'%.4f')
        "bulkTemperature=" + num2str(Tb,'%.2f')
        "wallTemperature=" + num2str(Tw,'%.2f')
        "fluid=Air"
        "viscosity=Sutherland"
        "dynamicInput=110.4"
    ], '&');

    %% ---- POST /generate_plot ---- %%

    post_headers = [
        HeaderField('Content-Type','application/x-www-form-urlencoded')
        HeaderField('Referer',[base_url '/input'])
        HeaderField('Origin',base_url)
    ].';

    post_req = RequestMessage( ...
        'POST', ...
        post_headers, ...
        char(body));

    post_resp = send(post_req, URI([base_url '/generate_plot']));

    if post_resp.StatusCode ~= 200
        error('thermoturb_query:postFailed', ...
              'POST failed with status %d', double(post_resp.StatusCode));
    end

    %% ---- Extract session cookie ---- %%

    set_cookie_field = post_resp.getFields('Set-Cookie');
    if isempty(set_cookie_field)
        error('thermoturb_query:noCookie', 'No Set-Cookie header found.');
    end

    cookie_str   = set_cookie_field.Value;
    cookie_value = regexp(cookie_str, '^[^;]+', 'match', 'once');

    %% ---- GET /download_data ---- %%

    cb = num2str(round(posixtime(datetime('now')) * 1000));

    get_headers = HeaderField('Cookie', cookie_value);

    get_req = RequestMessage( ...
        'GET', ...
        get_headers);

    get_resp = send(get_req, URI([base_url '/download_data?cb=' cb]));

    if get_resp.StatusCode ~= 200
        error('thermoturb_query:getFailed', ...
              'GET failed with status %d', double(get_resp.StatusCode));
    end

    %% ---- Unzip response ---- %%

    zip_bytes = get_resp.Body.Data;

    tmp_zip = [tempname() '.zip'];
    fid = fopen(tmp_zip, 'wb');
    fwrite(fid, zip_bytes);
    fclose(fid);

    tmp_dir = tempname();
    mkdir(tmp_dir);
    unzip(tmp_zip, tmp_dir);

    %% ---- Parse coefficients.dat ---- %%

    dat_path = fullfile(tmp_dir, 'coefficients.dat');
    raw = readlines(dat_path);
    raw = raw(strlength(strtrim(raw)) > 0);

    last_data_line = raw(end);
    vals = sscanf(last_data_line, '%f');

    % Columns: Retau, Rebulk, Cf, Stanton, Nusselt, Retau_cp
    Re_b   = vals(2);
    Cf_out = vals(3);
    Nu_out = vals(5);

    %% ---- Cleanup ---- %%

    delete(tmp_zip);
    rmdir(tmp_dir, 's');

end