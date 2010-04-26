/** Stereographic projection class. */
function StereographicProjection (phi1, lam1, rad) {
    this.phi1 = phi1;
    this.lam1 = lam1;
    this.rad = rad;
}

StereographicProjection.prototype.setCoords = function (phi1, lam1) {
    this.phi1 = phi1;
    this.lam1 = lam1;
}

StereographicProjection.prototype.setRadius = function (rad) {
    this.rad = rad;
}

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
    var cphi = Math.cos(phi1), sphi = Math.sin(phi1);
    var clam = Math.cos(lam1), slam = Math.sin(lam1);
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

    var cp1 = Math.cos(phi1);
    var dlam = lam-lam1;
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
            'r': rho
        };
    }
}

StereographicProjection.prototype.projectParallel = function (phi) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var R = this.rad;

    var s = Math.sin(phi1) + Math.sin(phi);
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
            'y': R*Math.cos(phi1)/s,
            'r': R*Math.cos(phi)/s
        };
    }
}

StereographicProjection.prototype.projectObj = function (re, de) {
    var lam1 = this.lam1;
    var phi1 = this.phi1;
    var rad = this.rad;
    var cphi = Math.cos(phi1), sphi = Math.sin(phi1);
    de = lam1-de;
    var cosc = Math.cos(re), sinc = Math.sin(re);
    var cosl = Math.cos(de), sinl = Math.sin(de);
    var k = rad / (1.0 + sphi * sinc + cphi * cosc * cosl);
    var x = k * cosc * sinl, y = k * (cphi * sinc - sphi * cosc * cosl);
    return [x, y, x*x + y*y < rad*rad];
}

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
}

StarMap.Moon = function (size, color) {
    this.size = size;
    this.color = color;
}

StarMap.Moon.prototype.pl = { name: 'Moon' };

StarMap.Moon.prototype.getCoord = function (jct, earthPos, equ2ecl) {
    // earthPos and equ2ecl are ignored
    var pos = StarJs.Solar.approxMoon(jct);
    return {'phi': pos.ra, 'theta': pos.dec};
}

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
    for (i = -180; i < 180; i += 20) {
        var p = this.proj.projectMeridian(Math.PI*i/180);
        ctx.beginPath();
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

    // Constellations
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(255,255,255,0.4)';
    for (j = clen; j--; ) {
        var s = co[j][0], e = co[j][1];
        var so = ortho[s], eo = ortho[e];
        if (so[3] || eo[3]) {
            ctx.moveTo((so[1]+halfsize), (halfsize-so[2]));
            ctx.lineTo((eo[1]+halfsize), (halfsize-eo[2]));
        }
    }
    ctx.stroke();

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
                var xx = cm[0]+halfsize, yy = halfsize-cm[1];
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
                var xx = cm[0]+halfsize, yy = halfsize-cm[1];
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
};

window['StarMap']=StarMap;
StarMap.prototype['setPos'] = StarMap.prototype.setPos;
