var LKS92WGS84 = (function()
{
    // Koordinātu pārveidojumos izmantotās konstantes
    LKS92WGS84.PI = Math.PI;                                    // Skaitlis pi
    LKS92WGS84.A_AXIS = 6378137;                                // Elipses modeļa lielā ass (a)
    LKS92WGS84.B_AXIS = 6356752.31414;                          // Elipses modeļa mazā ass (b)
    LKS92WGS84.CENTRAL_MERIDIAN = LKS92WGS84.PI * 24 / 180;     // Centrālais meridiāns
    LKS92WGS84.OFFSET_X = 500000;                               // Koordinātu nobīde horizontālās (x) ass virzienā
    LKS92WGS84.OFFSET_Y = -6000000;                             // Koordinātu nobīde vertikālās (y) ass virzienā
    LKS92WGS84.SCALE = 0.9996;                                  // Kartes mērogojuma faktors (reizinātājs)

    function LKS92WGS84() {}

    // Aprēķina loka garumu no ekvatora līdz dotā punkta ģeogrāfiskajam platumam
    LKS92WGS84.getArcLengthOfMeridian = function(phi)
    {
        var alpha, beta, gamma, delta, epsilon, n;

        n = (LKS92WGS84.A_AXIS - LKS92WGS84.B_AXIS) / (LKS92WGS84.A_AXIS + LKS92WGS84.B_AXIS);
        alpha = ((LKS92WGS84.A_AXIS + LKS92WGS84.B_AXIS) / 2) * (1 + (Math.pow(n, 2) / 4) + (Math.pow(n, 4) / 64));
        beta = (-3 * n / 2) + (9 * Math.pow(n, 3) / 16) + (-3 * Math.pow(n, 5) / 32);
        gamma = (15 * Math.pow(n, 2) / 16) + (-15 * Math.pow(n, 4) / 32);
        delta = (-35 * Math.pow(n, 3) / 48) + (105 * Math.pow(n, 5) / 256);
        epsilon = (315 * Math.pow(n, 4) / 512);

        return alpha * (phi + (beta * Math.sin(2 * phi)) + (gamma * Math.sin(4 * phi)) + (delta * Math.sin(6 * phi)) + (epsilon * Math.sin(8 * phi)));
    };

    // Aprēķina ģeogrāfisko platumu centrālā meridiāna punktam
    LKS92WGS84.getFootpointLatitude = function(y)
    {
        var yd, alpha, beta, gamma, delta, epsilon, n;

        n = (LKS92WGS84.A_AXIS - LKS92WGS84.B_AXIS) / (LKS92WGS84.A_AXIS + LKS92WGS84.B_AXIS);
        alpha = ((LKS92WGS84.A_AXIS + LKS92WGS84.B_AXIS) / 2) * (1 + (Math.pow(n, 2) / 4) + (Math.pow(n, 4) / 64));
        yd = y / alpha;
        beta = (3 * n / 2) + (-27 * Math.pow(n, 3) / 32) + (269 * Math.pow(n, 5) / 512);
        gamma = (21 * Math.pow(n, 2) / 16) + (-55 * Math.pow(n, 4) / 32);
        delta = (151 * Math.pow(n, 3) / 96) + (-417 * Math.pow(n, 5) / 128);
        epsilon = (1097 * Math.pow(n, 4) / 512);

        return yd + (beta * Math.sin(2 * yd)) + (gamma * Math.sin(4 * yd)) + (delta * Math.sin(6 * yd)) + (epsilon * Math.sin(8 * yd));
    };

    // Pārveido punkta ģeogrāfiskā platuma, garuma koordinātas par x, y koordinātām (bez pārvietojuma un mērogojuma)
    LKS92WGS84.convertMapLatLngToXY = function(phi, lambda, lambda0)
    {
        var N, nu2, ep2, t, t2, l,
            l3coef, l4coef, l5coef, l6coef, l7coef, l8coef,
            xy = [0, 0];

        ep2 = (Math.pow(LKS92WGS84.A_AXIS, 2) - Math.pow(LKS92WGS84.B_AXIS, 2)) / Math.pow(LKS92WGS84.B_AXIS, 2);
        nu2 = ep2 * Math.pow(Math.cos(phi), 2);
        N = Math.pow(LKS92WGS84.A_AXIS, 2) / (LKS92WGS84.B_AXIS * Math.sqrt(1 + nu2));
        t = Math.tan(phi);
        t2 = t * t;

        l = lambda - lambda0;
        l3coef = 1 - t2 + nu2;
        l4coef = 5 - t2 + 9 * nu2 + 4 * (nu2 * nu2);
        l5coef = 5 - 18 * t2 + (t2 * t2) + 14 * nu2 - 58 * t2 * nu2;
        l6coef = 61 - 58 * t2 + (t2 * t2) + 270 * nu2 - 330 * t2 * nu2;
        l7coef = 61 - 479 * t2 + 179 * (t2 * t2) - (t2 * t2 * t2);
        l8coef = 1385 - 3111 * t2 + 543 * (t2 * t2) - (t2 * t2 * t2);

        // x koordināta
        xy[0] = N * Math.cos(phi) * l + (N / 6 * Math.pow(Math.cos(phi), 3) * l3coef * Math.pow(l, 3)) + (N / 120 * Math.pow(Math.cos(phi), 5) * l5coef * Math.pow(l, 5)) + (N / 5040 * Math.pow(Math.cos(phi), 7) * l7coef * Math.pow(l, 7));

        // y koordināta
        xy[1] = LKS92WGS84.getArcLengthOfMeridian(phi) + (t / 2 * N * Math.pow(Math.cos(phi), 2) * Math.pow(l, 2)) + (t / 24 * N * Math.pow(Math.cos(phi), 4) * l4coef * Math.pow(l, 4)) + (t / 720 * N * Math.pow(Math.cos(phi), 6) * l6coef * Math.pow(l, 6)) + (t / 40320 * N * Math.pow(Math.cos(phi), 8) * l8coef * Math.pow(l, 8));

        return xy;
    };

    // Pārveido punkta x, y koordinātas par ģeogrāfiskā platuma, garuma koordinātām (bez pārvietojuma un mērogojuma)
    LKS92WGS84.convertMapXYToLatLon = function(x, y, lambda0)
    {
        var phif, Nf, Nfpow, nuf2, ep2, tf, tf2, tf4, cf,
            x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac,
            x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly,
            latLng = [0, 0];

        phif = LKS92WGS84.getFootpointLatitude(y);
        ep2 = (Math.pow(LKS92WGS84.A_AXIS, 2) - Math.pow(LKS92WGS84.B_AXIS, 2)) / Math.pow(LKS92WGS84.B_AXIS, 2);
        cf = Math.cos(phif);
        nuf2 = ep2 * Math.pow(cf, 2);
        Nf = Math.pow(LKS92WGS84.A_AXIS, 2) / (LKS92WGS84.B_AXIS * Math.sqrt(1 + nuf2));
        Nfpow = Nf;

        tf = Math.tan(phif);
        tf2 = tf * tf;
        tf4 = tf2 * tf2;

        x1frac = 1 / (Nfpow * cf);

        Nfpow *= Nf;    // Nf^2
        x2frac = tf / (2 * Nfpow);

        Nfpow *= Nf;    // Nf^3
        x3frac = 1 / (6 * Nfpow * cf);

        Nfpow *= Nf;    // Nf^4
        x4frac = tf / (24 * Nfpow);

        Nfpow *= Nf;    // Nf^5
        x5frac = 1 / (120 * Nfpow * cf);

        Nfpow *= Nf;    // Nf^6
        x6frac = tf / (720 * Nfpow);

        Nfpow *= Nf;    // Nf^7
        x7frac = 1 / (5040 * Nfpow * cf);

        Nfpow *= Nf;    // Nf^8
        x8frac = tf / (40320 * Nfpow);

        x2poly = -1 - nuf2;
        x3poly = -1 - 2 * tf2 - nuf2;
        x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 - 3 * (nuf2 * nuf2) - 9 * tf2 * (nuf2 * nuf2);
        x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2;
        x6poly = -61 - 58 * tf2 + (tf2 * tf2) + 270 * nuf2 - 330 * tf2 * nuf2;
        x7poly = 61 - 479 * tf2 + 179 * (tf2 * tf2) - (tf2 * tf2 * tf2);
        x8poly = 1385 - 3111 * tf2 + 543 * (tf2 * tf2) - (tf2 * tf2 * tf2);

        latLng[0] = phif + x2frac * x2poly * (x * x) + x3frac * x3poly * (x * x * x) + x4frac * x4poly * (x * x * x * x) + x5frac * x5poly * (x * x * x * x * x) + x6frac * x6poly * (x * x * x * x * x * x) + x7frac * x7poly * (x * x * x * x * x * x * x) + x8frac * x8poly * (x * x * x * x * x * x * x * x);
        latLng[1] = lambda0 + x1frac * x + x2frac * x2poly * (x * x) + x3frac * x3poly * (x * x * x) + x4frac * x4poly * (x * x * x * x) + x5frac * x5poly * (x * x * x * x * x) + x6frac * x6poly * (x * x * x * x * x * x) + x7frac * x7poly * (x * x * x * x * x * x * x) + x8frac * x8poly * (x * x * x * x * x * x * x * x);

        return latLng;
    };

    // Pārveido punkta ģeogrāfiskā platuma, garuma koordinātas par x, y koordinātām (ar pārvietojumu un mērogojumu)
    LKS92WGS84.convertLatLonToXY = function(coordinates)
    {
        var phi = coordinates[0] * LKS92WGS84.PI / 180;
        var lambda = coordinates[1] * LKS92WGS84.PI / 180;
        var xy = LKS92WGS84.convertMapLatLngToXY(phi, lambda, LKS92WGS84.CENTRAL_MERIDIAN);
        xy[0] = xy[0] * LKS92WGS84.SCALE + LKS92WGS84.OFFSET_X;
        xy[1] = xy[1] * LKS92WGS84.SCALE + LKS92WGS84.OFFSET_Y;
        return xy;
    };

    // Pārveido punkta x, y koordinātas par ģeogrāfiskā platuma, garuma koordinātām (ar pārvietojumu un mērogojumu)
    LKS92WGS84.convertXYToLatLon = function(coordinates)
    {
        var x = (coordinates[0] - LKS92WGS84.OFFSET_X) / LKS92WGS84.SCALE;
        var y = (coordinates[1] - LKS92WGS84.OFFSET_Y) / LKS92WGS84.SCALE;
        var latLng = LKS92WGS84.convertMapXYToLatLon(x, y, LKS92WGS84.CENTRAL_MERIDIAN);
        latLng[0] = latLng[0] * 180 / LKS92WGS84.PI;
        latLng[1] = latLng[1] * 180 / LKS92WGS84.PI;
        return latLng;
    };

    return LKS92WGS84;
})();
// Define the EPSG:3059 projection
proj4.defs("EPSG:3059", "+proj=tmerc +lat_0=0 +lon_0=24 +k=0.9996 +x_0=500000 +y_0=-6000000 +datum=WGS84 +units=m +no_defs");

// Initialize the map
var map = L.map('map').setView([56.8796, 24.6032], 8); // Centered on Latvia

// Add the OpenStreetMap tiles
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
}).addTo(map);

// Function to convert coordinates from EPSG:3059 to EPSG:4326
function convertCoordinates(coords) {
    return proj4('EPSG:3059', 'EPSG:4326', coords);
}

// Load GeoJSON data
fetch('file.json')
    .then(response => response.json())
    .then(data => {
        // Convert coordinates in the GeoJSON data
        data.features.forEach(feature => {
            if (feature.geometry.type === 'Point') {
                feature.geometry.coordinates = convertCoordinates(feature.geometry.coordinates);
                console.log('Converted coordinates:', feature.geometry.coordinates); // Log converted coordinates
            }
        });

        // Add GeoJSON layer to the map
        L.geoJSON(data, {
            onEachFeature: function (feature, layer) {
                if (feature.properties) {
                    var popupContent = '<strong>' + feature.properties.PLACENAME + '</strong><br>';
                    for (var property in feature.properties) {
                        if (feature.properties.hasOwnProperty(property) && property !== 'PLACENAME') {
                            popupContent += property + ': ' + feature.properties[property] + '<br>';
                        }
                    }
                    layer.bindPopup(popupContent);
                }
            },
            pointToLayer: function (feature, latlng) {
                return L.marker(latlng);
            }
        }).addTo(map);
    })
    .catch(error => console.error('Error loading the GeoJSON data:', error));