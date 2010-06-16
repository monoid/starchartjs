/** Stereographic projection class. */
function StereographicProjection (phi1, lam1, rad) {
    this.setCoords(phi1, lam1);
    this.setRadius(rad);
}

StereographicProjection.prototype.setCoords = function (phi1, lam1) {
    this.phi1 = phi1;
    this.lam1 = lam1;

    this.cph1 = Math.cos(phi1);
    this.sph1 = Math.sin(phi1);
    
    this.cla1 = Math.cos(lam1);
    this.sla1 = Math.sin(lam1);
    
};

StereographicProjection.prototype.setRadius = function (rad) {
    this.rad = rad;
};

StereographicProjection.prototype.projectPoints = function (arr, rad) {
    function sinSum(cosa, sina, cosb, sinb) {
        return cosa*sinb+sina*cosb;
    }

    function cosSum(cosa, sina, cosb, sinb) {
        return cosa*cosb-sina*sinb;
    }

    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var rad = this.rad;

    var len = arr.length, i;
    var res = Array(len);
    var cphi = this.cph1, sphi = this.sph1;
    var clam = this.cla1, slam = this.sla1;
    for (i = 0; i < len; ++i) {
        var star = arr[i];
        var mag = star[0], re = star[2], de = -star[1];
        var t2c = re*re, t2l = de*de, t2c1=1+t2c, t2l1=1+t2l;
        var cosc = (1-t2c)/t2c1, sinc = 2*re/t2c1;
        var cosl1 = (1-t2l)/t2l1, sinl1 = 2*de/t2l1;
        var cosl = cosSum(cosl1, sinl1, clam, slam), sinl = sinSum(cosl1, sinl1, clam, slam);
        var k = rad / (1.0 + sphi * sinc + cphi * cosc * cosl);
        var x = k * cosc * sinl, y = k * (cphi * sinc - sphi * cosc * cosl);
        res[i] = [mag,
                  x,
                  y,
                  x*x + y*y < rad*rad];
    }
    return res;
};

StereographicProjection.prototype.projectMeridian = function (lam) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var R = this.rad;

    var cp1 = this.cph1;
    var dlam = lam1-lam;
    var sl1 = Math.sin(dlam);
    if (Math.abs(sl1) < 1e-10 || Math.abs(cp1) < 1e-10) {
        return {
            'type': 'line',
            'x': 0,
            'y': 0,
            'vx': Math.cos(dlam),
            'vy': sl1
        };
    } else {
        var x = -R/(cp1*Math.tan(dlam));
        var y = -R*Math.tan(phi1);
        var rho = R/(cp1*sl1);
        return {
            'type': 'circle',
            'x': x,
            'y': y,
            'r': Math.abs(rho)
        };
    }
};

StereographicProjection.prototype.projectParallel = function (phi) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var R = this.rad;

    var s = this.sph1 + Math.sin(phi);
    // TODO: line if s == 0
    if (Math.abs(s) < 1e-10) {
        return {
            'type': 'line',
            'x': 0,
            'y': 0,
            'vx': 1,
            'vy': 0
        };
    } else {
        return {
            'type': 'circle',
            'x': 0,
            'y': R*this.cph1/s,
            'r': Math.abs(R*Math.cos(phi)/s)
        };
    }
};

StereographicProjection.prototype.projectObj = function (re, de) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var rad = this.rad;
    var cphi = this.cph1, sphi = this.sph1;
    de = lam1-de;
    var cosc = Math.cos(re), sinc = Math.sin(re);
    var cosl = Math.cos(de), sinl = Math.sin(de);
    var k = rad / (1.0 + sphi * sinc + cphi * cosc * cosl);
    var x = k * cosc * sinl, y = - k * (cphi * sinc - sphi * cosc * cosl);
    return [x, y, x*x + y*y < rad*rad];
};

StereographicProjection.prototype.projectSegment = function (ra1, de1, ra2, de2) {
/*

Produced by Maxima 5.18.1:

(%i1) z1: re1 + %i*im1;
(%o1)                            re1 + %i im1
(%i2) z2: re2 + %i*im2;
(%o2)                            re2 + %i im2
(%i3) z0: ((abs(z1)^2-1)*z2 - (abs(z2)^2-1)*z1)/(conjugate(z1)*z2 - z1*conjugate(z2));
          2      2                                          2      2
      (re1  + im1  - 1) (re2 + %i im2) - (re1 + %i im1) (re2  + im2  - 1)
(%o3) -------------------------------------------------------------------
         (re1 - %i im1) (re2 + %i im2) - (re1 + %i im1) (re2 - %i im2)
(%i4) realpart(z0);
                         2      2                2      2
                 im2 (re1  + im1  - 1) - im1 (re2  + im2  - 1)
(%o4)            ---------------------------------------------
                             2 im2 re1 - 2 im1 re2
(%i5) imagpart(z0);
                      2      2                    2      2
                  (re1  + im1  - 1) re2 - re1 (re2  + im2  - 1)
(%o5)           - ---------------------------------------------
                              2 im2 re1 - 2 im1 re2
*/

    function s(x) { return x*x; }

    var p1 = this.projectObj(ra1, de1), p2 = this.projectObj(ra2, de2);

    var re1 = p1[0]/this.rad, im1 = p1[1]/this.rad;
    var re2 = p2[0]/this.rad, im2 = p2[1]/this.rad;

    var a1 = re1*re1 + im1*im1 - 1;
    var a2 = re2*re2 + im2*im2 - 1;

    var ha1 =   a1 * im2 - a2 * im1;
    var ha2 = - a1 * re2 + a2 * re1;

    var d = 2*(im2 * re1 - im1 * re2);
    // find center
    var c1 = ha1 / d;
    var c2 = ha2 / d;

    // find radius
    var r = Math.sqrt(1 + c1*c1 + c2*c2);
    // find angles
    var a1 = Math.atan2(im1 - c2, re1 - c1);
    var a2 = Math.atan2(im2 - c2, re2 - c1);
    return {
        'type': 'circle',
        'x': c1*this.rad, 'y': c2*this.rad,
        'a1': a1, 'a2': a2,
        'p1': p1, 'p2': p2,
        'r': r*this.rad
    };
};


StereographicProjection.prototype.projectParallelSegment = function (lam1, lam2, phi) {
    var p = this.projectParallel(phi);
    var p1 = this.projectObj(phi, lam1), p2 = this.projectObj(phi, lam2);
    if (p.type === 'circle') {
        p.y = - p.y;
        p.a1 = Math.atan2(p1[1] - p.y, p1[0] - p.x);
        p.a2 = Math.atan2(p2[1] - p.y, p2[0] - p.x);
        return p;
    } else {
        return {
            type: 'line',
            x1: p1[0], y1: p1[1],
            x2: p2[0], y2: p2[1]
        };
    }
};


/** Project circle centered at point p (already projected to plane)
 * with angular radius alpha.
 *
 * p is usually result of StereographicProjection.prototype.projectObj.
 */
StereographicProjection.prototype.projectCircle2 = function (p, alpha) {
    var r = Math.tan(alpha/2);

    var aa2 = (p[0]*p[0]+p[1]*p[1])/(this.rad*this.rad);
    var denom = 1 - aa2*r*r;

    var rad = r*(aa2+1)/denom;
    var cx = p[0]*(1+r*r)/denom;
    var cy = p[1]*(1+r*r)/denom;

    return {
        'type': 'circle',
        'x': cx, 'y': cy,
        'rad': this.rad*rad
    };
};

/** Project circle centered at (re, de) on the sphere with angular
 * radius alpha.
 */
StereographicProjection.prototype.projectCircle = function (re, de, alpha) {
    var p = this.projectObj(re, de);
    return this.projectCircle2(p, alpha);
};

/** Inverse projection of a single point.
 */
StereographicProjection.prototype.inverseObj = function (x, y) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;

    var cp1 = this.cph1, sp1 = this.sph1;

    var rho2 = x*x + y*y;
    if (rho2 > this.rad*this.rad) {
        // Point is out of map
        return null;
    }
    if (rho2 < 1e-18) {
        return {
            're': lam1,
            'de': phi1
        };
    } else {
        var rho = Math.sqrt(rho2);
        var c = 2*Math.atan(rho/this.rad);
        var cosc = Math.cos(c);
        var sinc = Math.sin(c);
        var phi = Math.asin(cosc*sp1
                            + (y*sinc*cp1/rho));
        var lam;
        if (Math.abs(phi1-Math.PI) < 1e-9) {
            // phi1 = PI
            lam = lam1 - Math.atan2(x, -y);
        } else if (Math.abs(this.phi - (-Math.PI)) < 1e-9) {
            // phi1 = -PI
            lam = lam1 - Math.atan2(x, y);
        } else {
            lam = lam1 - Math.atan2(x*sinc,
                                    rho*cp1*cosc - y*sp1*sinc);
        }
        lam = StarJs.Math.mod(lam, 2*Math.PI);
        phi = StarJs.Math.mod(phi+Math.PI, 2*Math.PI)-Math.PI;
        return {
            're': lam,
            'de': phi
        };
    }
};

/**
 * Celestial map component.
 * @constructor
 */
function StarMap (elt, size, stars, cnstltns, prop) {
    this.paper = document.getElementById(elt);
    this.ctx = this.paper.getContext("2d");
    this.prop = prop;

    this.size = size;
    var halfsize = Math.floor(size/2);

    this.planets = (typeof prop.planets === 'undefined') ? true : prop.planets;

    this.stars = stars;
    this.cnstltns = cnstltns;

    this.proj = new StereographicProjection(0, 0, halfsize);

    this.drawBg();
}

StarMap.prototype.drawBg = function () {
    var size = this.size;
    var halfsize = Math.floor(size/2);
    var ctx = this.ctx;

    ctx.clearRect(0, 0, size, size);
    ctx.fillStyle='#FFF';
    ctx.fillRect(0,0,size, size);
    ctx.beginPath();
    ctx.fillStyle = (this.prop.circleFill || "#000010");
    ctx.arc(halfsize, halfsize, halfsize, 0, 2*Math.PI, true);
    ctx.fill();
    
    ctx.beginPath();
    ctx.arc(halfsize, halfsize, halfsize, 0, 2*Math.PI, true);
    ctx.clip();
}    

StarMap.Planet = function (pl, size, color) {
    this.pl = pl;
    this.size = size;
    this.color = color;
}

StarMap.Planet.prototype.getCoord = function (jct, earthPos, equ2ecl) {
    var pos = this.pl.keplerCoord(jct);
    return new StarJs.Vector.Polar3(equ2ecl.apply(pos.sub(earthPos)))
};

StarMap.Moon = function (size, color) {
    this.size = size;
    this.color = color;
};

StarMap.Moon.prototype.pl = { name: 'Moon' };

StarMap.Moon.prototype.getCoord = function (jct, earthPos, equ2ecl) {
    // earthPos and equ2ecl are ignored
    var pos = StarJs.Solar.approxMoon(jct);
    return {'phi': pos.ra, 'theta': pos.dec};
};

StarMap.PLANETS = [
    new StarMap.Planet(StarJs.Solar.BODIES.Sun, 20, '#FF0'),
    new StarMap.Moon(20, '#880'),
    new StarMap.Planet(StarJs.Solar.BODIES.Mercury, 3, '#888'),
    new StarMap.Planet(StarJs.Solar.BODIES.Venus, 4, '#AAA'),
    new StarMap.Planet(StarJs.Solar.BODIES.Mars, 4, '#F80'),
    new StarMap.Planet(StarJs.Solar.BODIES.Jupiter, 6, '#FB0'),
    new StarMap.Planet(StarJs.Solar.BODIES.Saturn, 6, '#AA0'),
    new StarMap.Planet(StarJs.Solar.BODIES.Uranus, 6, '#CAF'),
    new StarMap.Planet(StarJs.Solar.BODIES.Neptune, 6, '#CAF')
];

StarMap.EARTH = StarJs.Solar.BODIES.Earth;

StarMap.prototype.drawTelrad = function (ctx, lat, lon) {
    var p = this.proj.projectObj(lat, lon);
    var g05 = this.proj.projectCircle2(p, 0.5/180*Math.PI);
    var g20 = this.proj.projectCircle2(p, 2.0/180*Math.PI);
    var g40 = this.proj.projectCircle2(p, 4.0/180*Math.PI);

    var h = Math.floor(this.size/2);

    ctx.strokeStyle = 'rgba(255,0,0,0.6)';
    ctx.lineWidth = g40.rad/9;

    function drawBullEye(g) {
        var D = 0.15;
        ctx.beginPath();
        ctx.arc(h+g.x, h+g.y, g.rad, D, 0.5*Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(h+g.x, h+g.y, g.rad, 0.5*Math.PI+D, Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(h+g.x, h+g.y, g.rad, Math.PI+D, 1.5*Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(h+g.x, h+g.y, g.rad, 1.5*Math.PI+D, 2*Math.PI-D, false);
        ctx.stroke();
    }

    ctx.beginPath();
    ctx.arc(h+g05.x, h+g05.y, g05.rad, 0, 2*Math.PI, true);
    ctx.stroke();

    drawBullEye(g20);
    drawBullEye(g40);
};

StarMap.prototype.setPos = function (lat, lon, time) {
    var Ti = StarJs.Time;

    if (typeof time === 'undefined') {
        time = +new Date();
    } else if (typeof time !== 'number') {
        time = +time;
    }
    
    var mjd = Ti.time2mjd(time);
    var gms_t = Ti.gmst(mjd);

    /** @const */
    var DEG2RAD = StarJs.Math.DEG2RAD;
    lat *= DEG2RAD;
    lon *= DEG2RAD;

    lat += gms_t;

    this.proj.setCoords(lon, lat);

    var ortho = this.proj.projectPoints(this.stars);
    var cst = [], i, j, slen = ortho.length, co = this.cnstltns, clen = co.length, halfsize = Math.floor(this.size/2);
    
    this.drawBg();

    var ctx = this.ctx;

    // Draw graticule
    ctx.strokeStyle = '#448';
    for (i = -80; i < 90; i += 10) {
        var p = this.proj.projectParallel(Math.PI*i/180);
        ctx.beginPath();
        ctx.lineWidth = (i === 0) ? 1.7 : 1;
        switch (p.type) {
        case 'line':
            ctx.moveTo(halfsize+p.x-halfsize*p.vx,
                       halfsize+p.y-halfsize*p.vy);
            ctx.lineTo(halfsize+p.x+halfsize*p.vx,
                       halfsize+p.y+halfsize*p.vy);
            break;
        case 'circle':
            ctx.arc(halfsize+p.x, halfsize-p.y, p.r, 0, 2*Math.PI, true);
            break;
        }
        ctx.stroke();
    }
    for (i = 0; i < 180; i += 15) {
        var p = this.proj.projectMeridian(Math.PI*i/180);
        ctx.beginPath();
        ctx.lineWidth = (i === 0 || i === -180) ? 1.7 : 1;
        switch (p.type) {
        case 'line':
            ctx.moveTo(halfsize+p.x-halfsize*p.vx,
                       halfsize+p.y-halfsize*p.vy);
            ctx.lineTo(halfsize+p.x+halfsize*p.vx,
                       halfsize+p.y+halfsize*p.vy);
            break;
        case 'circle':
            ctx.arc(halfsize+p.x, halfsize-p.y, p.r,
                    0, 2*Math.PI, true);
            break;
        }
        ctx.stroke();
    }
    ctx.lineWidth = 1;

    function angSep(a1, a2) {
        return StarJs.Math.mod(a1-a2, 2 * Math.PI);
    };
    // Constellation boundaries
    var cstn = null;
    // TODO: make a = 0.5 after boundaries line merging.  Now many
    // lines are drawn twice, and they color differs (somewhat
    // brighter) if alpha is used.
    ctx.strokeStyle = 'rgba(128,0,128,1)';
    var prev = null;
    for (j = 0; j < CON_BOUND_18.length; ++j) {
        var l = CON_BOUND_18[j];
        if (cstn === l[2]) {
            var seg;
            if (prev[1] === l[1]) {
                seg = this.proj.projectParallelSegment(
                    15*Math.PI*prev[0]/180,
                    15*Math.PI*l[0]/180,
                    Math.PI*l[1]/180
                );
            } else {
                seg = this.proj.projectSegment(Math.PI*prev[1]/180,
                                               15*Math.PI*prev[0]/180,
                                               Math.PI*l[1]/180,
                                               15*Math.PI*l[0]/180);
            }
            ctx.beginPath();
            switch (seg.type) {
            case 'circle':   
                ctx.arc(halfsize+seg.x, halfsize+seg.y,
                        seg.r,
                        seg.a1, seg.a2,
                        // TODO: angular sepration of unprojected
                        // lines.  Projected circles sometimes give
                        // wrong result due to distortion.
                        angSep(seg.a1,seg.a2) < angSep(seg.a2,seg.a1));
                break;
            case 'line':
                // TODO sometimes lines shouldn't be drawn if their
                // central point pass through infinity, or at least
                // they should be drawn more intelligently.  Can be
                // such lines wisible in viewport?
                ctx.moveTo(halfsize+seg.x1, halfsize+seg.y1);
                ctx.lineTo(halfsize+seg.x2, halfsize+seg.y2);
                break;
            }
            ctx.stroke();
        } else {
            // TODO: draw a closing arc.
            prev = l;
            cstn = l[2];
        }
        prev = l;
    }

    // Constellations
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(255,255,255,0.7)';
    for (j = clen; j--; ) {
        var s = co[j][0], e = co[j][1];
        var so = ortho[s], eo = ortho[e];
        if (so[3] || eo[3]) {
            ctx.moveTo((so[1]+halfsize), (halfsize-so[2]));
            ctx.lineTo((eo[1]+halfsize), (halfsize-eo[2]));
        }
    }
    ctx.stroke();

    // Draw ecliptics
    var eclp = this.proj.projectCircle(Math.PI/2 + StarJs.Solar.EPS, Math.PI/2, Math.PI/2);
    if (eclp.type === 'circle') {
        ctx.beginPath();
        ctx.strokeStyle = 'yellow';
        ctx.arc(eclp.x+halfsize, halfsize+eclp.y, eclp.rad, 0, 2*Math.PI, true);
        ctx.stroke();
    }

    // Stars
    ctx.fillStyle = '#FFF';
    for (i = 0; i < slen; ++i) {
        var s = ortho[i];
        if (s[3]) {
            ctx.beginPath();
            ctx.arc(s[1]+halfsize, halfsize-s[2],
                    Math.max(3-s[0]/2, 0.5),
                    0, 2*Math.PI, true);
            ctx.fill();
        }
    }

    function messierColor(mag) {
        var v = Math.min(15, Math.floor(19-mag));
        var h = v.toString(16);
        return '#'+h+h+h;
    }

    // Draw Messier objects
    if (this.prop && this.prop.messier) {
        var messier = this.prop.messier;
        var mlen = messier.length, cc, cm;
        for (i = 0; i < mlen; ++i) {
            cc = messier[i];
            cm = this.proj.projectObj(cc[4], 15*cc[3]);
            if (cm[2]) {
                var xx = cm[0]+halfsize, yy = halfsize+cm[1];
                ctx.beginPath();
                if (this.prop.messier_colors && this.prop.messier_colors[messier[i][2]]) {
                    ctx.strokeStyle = this.prop.messier_colors[cc[2]];
                } else {
                    ctx.strokeStyle = messierColor(cc[5]);
                }
                ctx.arc(xx, yy, 4, 0, 2*Math.PI, true);
                ctx.stroke();
            }
        }
    }

    // Draw planets
    if (this.planets) {
        var jct = Ti.mjd2jct(mjd);
        var earthPos = StarMap.EARTH.keplerCoord(jct);
        var equ2ecl = StarJs.Coord.ecl2equMatrix(jct);
        for (i = 0; i < StarMap.PLANETS.length; ++i) {
            var planet = StarMap.PLANETS[i];
            cc = planet.getCoord(jct, earthPos, equ2ecl);
            
            cm = this.proj.projectObj(cc.theta, cc.phi);
            if (cm[2]) {
                ctx.beginPath();
                ctx.fillStyle = planet.color;
                var xx = cm[0]+halfsize, yy = halfsize + cm[1];
                ctx.arc(xx, yy, planet.size/2,
                        0, 2*Math.PI, true);
                ctx.fill();
                ctx.beginPath();
                ctx.strokeStyle = planet.color;
                ctx.arc(xx, yy, planet.size/2 + 2,
                        0, 2*Math.PI, true);
                ctx.stroke();
                
                if (ctx.fillText) {
                    ctx.fillText(planet.pl.name,
                                 xx + planet.size/2 + 1,
                                 yy - planet.size/2 - 1);
                }
            }
        }
    }

    // Draw sample telrads
    var tel = {'x': 0.901, 'y': 0.451};
    this.drawTelrad(ctx, tel.x, tel.y);
    tel = {'x': 0.301, 'y': 2.151};
    this.drawTelrad(ctx, tel.x, tel.y);

    // Draw path of C/2009 R1 (McNaught)
    var C2009R1_DA = ["2010 05 04", "2010 05 09", "2010 05 14",
                      "2010 05 19", "2010 05 24", "2010 05 29",
                      "2010 06 03", "2010 06 08", "2010 06 13",
                      "2010 06 18", "2010 06 23", "2010 06 28",
                      "2010 07 03", "2010 07 08", "2010 07 13",
                      "2010 07 18", "2010 07 23", "2010 07 28"];
    var C2009R1_RA = [23+35.99/60.0, 23+48.09/60.0,  0+ 2.04/60.0,
                       0+18.56/60.0,  0+38.69/60.0,  1+ 4.06/60.0,
                       1+37.04/60.0,  2+20.73/60.0,  3+17.79/60.0,
                       4+26.49/60.0,  5+37.00/60.0,  6+36.84/60.0,
                       7+20.69/60.0,  7+50.73/60.0,  8+11.73/60.0,
                       8+27.57/60.0,  8+40.53/60.0,  8+51.84/60.0];
    var C2009R1_DE = [ 9+56.1/60.0, 13+11.8/60.0, 16+53.6/60.0,
                      21+ 5.1/60.0, 25+48.9/60.0, 31+ 3.2/60.0,
                      36+37.2/60.0, 42+ 1.6/60.0, 46+19.9/60.0,
                      48+14.0/60.0, 46+44.8/60.0, 42+ 0.7/60.0,
                      35+ 8.7/60.0, 27+30.3/60.0, 20+ 4.7/60.0,
                      13+16.4/60.0, 07+ 8.3/60.0,  1+35.8/60.0];
    ctx.fillStyle = 'green';
    ctx.strokeStyle = 'green';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (i = 0; i < C2009R1_DE.length; ++i) {
        var cp = this.proj.projectObj(Math.PI*C2009R1_DE[i]/180.0,
                                      Math.PI*C2009R1_RA[i]/12.0);
        if (cp[2]) {
            xx = cp[0]+halfsize;
            yy = halfsize-cp[1];
            ctx.lineTo(xx, yy);
        }
    }
    ctx.stroke();

    ctx.fillStyle = '#8F8';
    ctx.strokeStyle = '#8F8';
    // We recalculate same data, but this is just a sample, nevermind.
    for (i = 0; i < C2009R1_DE.length; ++i) {
        var cp = this.proj.projectObj(Math.PI*C2009R1_DE[i]/180.0,
                                      Math.PI*C2009R1_RA[i]/12.0);
        if (cp[2]) {
            xx = cp[0]+halfsize;
            yy = halfsize-cp[1];
            ctx.beginPath();
            ctx.arc(xx, yy, 2, 0, 2*Math.PI, true);
            ctx.fill();
            if (ctx.fillText) {
                ctx.fillText(C2009R1_DA[i], xx, yy - 4);
            }
        }
    }
};

window['StarMap']=StarMap;
StarMap.prototype['setPos'] = StarMap.prototype.setPos;
