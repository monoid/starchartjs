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

StereographicProjection.prototype.projectPoints = function (arr) {
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
        var x = k * cosc * sinl, y = k * (sphi * cosc * cosl - cphi * sinc);
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
            type: 'line',
            x: 0,
            y: 0,
            vx: Math.cos(dlam),
            vy: sl1
        };
    } else {
        var x = -R/(cp1*Math.tan(dlam));
        var y = R*Math.tan(phi1);
        var rho = R/(cp1*sl1);
        return {
            type: 'circle',
            x: x,
            y: y,
            flip: rho < 0,
	    r: Math.abs(rho)
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
            type: 'line',
            x: 0,
            y: 0,
            vx: 1,
            vy: 0
        };
    } else {
	var rho = R*Math.cos(phi)/s;
        return {
            type: 'circle',
            x: 0,
            y: -R*this.cph1/s,
            flip: rho < 0,
            r: Math.abs(rho)
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

/** Coordinates of great circle segment from (ra1,de1) to (ra2,de2).
 */
StereographicProjection.prototype.projectGreatSegment = function (ra1, de1, ra2, de2) {
/*

Produced by Maxima 5.18.1:

(%i1) z1: re1 + %i*im1;
(%o1)                            re1 + %i im1
(%i2) z2: re2 + %i*im2;
(%o2)                            re2 + %i im2

  Center of great circle defined by two points z1 and z2.  See (4.4) in

       Robin G. Stuart.  Application of Complex Analysis to Spherical
       Coordinate System.  Q.Jl. R. astr. Soc. (1984) 25, 126-136.

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
    a1 = Math.atan2(im1 - c2, re1 - c1);
    a2 = Math.atan2(im2 - c2, re2 - c1);
    return {
        type: 'circle',
        x: c1*this.rad, y: c2*this.rad,
        a1: a1, a2: a2,
        p1: p1, p2: p2,
        flip: d < 0,
        r: r*this.rad
    };
};


/** Project segment of parallel phi from lam1 to lam2.
 */
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

    var rho = r*(aa2+1)/denom;
    var cx = p[0]*(1+r*r)/denom;
    var cy = p[1]*(1+r*r)/denom;

    return {
        type: 'circle',
        x: cx, y: cy,
        flip: denom < 0,
        rad: this.rad*Math.abs(rho)
    };
};

/** Project circle centered at (re, de) on the sphere with angular
 * radius alpha.
 */
StereographicProjection.prototype.projectCircle = function (re, de, alpha) {
    var p = this.projectObj(re, de);
    return this.projectCircle2(p, alpha);
};


/** Project arc of arbitrary circle centered at (re, de) on the sphere
 * with angular radius alpha.  Arc spans from (re1, de1) to
 * (re2, de2).
 */

StereographicProjection.prototype.projectSegment = function (re, de, alpha, re1, de1, re2, de2) {
    var p1 = this.projectObj(re1, de1);
    var p2 = this.projectObj(re2, de2);

    var c = this.projectCircle(re, de, alpha);

    c.r = c.rad; // TODO Inconsistence...
    c.p1 = p1;
    c.p2 = p2;

    c.a1 = Math.atan2(p1[1]-c.y, p1[0]-c.x);
    c.a2 = Math.atan2(p2[1]-c.y, p2[0]-c.x);

    return c;
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
            'ra': lam1,
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
            'ra': lam,
            'de': phi
        };
    }
};
StereographicProjection.prototype['inverseObj'] = StereographicProjection.prototype.inverseObj;

/**
 * Rendering scene: objects are cached for rotation and zooming.  When
 * rotation or zoom is changed, internal matrix is recalculated, and
 * on render method this matrix is applied to all points.  We do not
 * use canvas' rotation and scaling because first affects labels (we
 * wat them horizontal) and second affects line with (we want
 * zoom-independent width).
 *
 * Possible workaround exist for this problem: keep orientation
 * internally and compensate it when rendering labels, and keep zoom
 * factor too and compensate it when drawing lines, arcs and circles.
 */
function Scene() {
    this.objects = [];
    this.matrix = null; // TODO
    this.scale = 1.0;
    this.rotation = 0.0;
}

Scene.prototype.drawLine = function (x1, y1, x2, y2, width, color) {
    this.objects.push(function (ctx) {
        ctx.beginPath();
        ctx.strokeColor = color;
        ctx.lineWidth = width/this.scale;
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
    });
};

Scene.prototype.drawArc = function (x, y, r, a1, a2, width, color) {
    this.objects.push(function (ctx) {
        var rot = this.rotation;
        ctx.beginPath();
        ctx.strokeColor = color;
        ctx.lineWidth = width/this.scale;
        ctx.arc(x, y, r, a1, a2);
        ctx.stroke();
    });
};

Scene.prototype.drawCircle = function (x, y, r, width, color) {
    this.objects.push(function (ctx) {
        var rot = this.rotation;
        ctx.beginPath();
        ctx.strokeColor = color;
        ctx.lineWidth = width/this.scale;
        ctx.arc(x, y, r, 0, 2*Math.PI);
        ctx.stroke();
    });
};

Scene.prototype.drawPlanet = function (x, y, r, title, color) {
};

Scene.prototype.drawStars = function (stars) {
};

Scene.prototype.reset = function () {
    this.objects = [];
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

    this.planets = (typeof prop['planets'] === 'undefined') ? true : prop['planets'];

    this.orient = 0;

    this.stars = stars;
    this.cnstltns = cnstltns;

    this['proj'] = new StereographicProjection(0, 0, halfsize);

    //this.drawBg();
}

StarMap.prototype.setSize = function (size) {
    this['proj'].setRadius(Math.floor(size/2));
    this.size = size;
};
StarMap.prototype['setSize'] = StarMap.prototype.setSize;

StarMap.prototype.drawBg = function () {
    var size = this.size;
    var halfsize = Math.floor(size/2);
    var ctx = this.ctx;

    ctx.clearRect(0, 0, size, size);
    ctx.fillStyle='#FFF';
    ctx.fillRect(-halfsize,-halfsize,size, size);

    ctx.beginPath();
    ctx.fillStyle = (this.prop['circleFill'] || "#000010");
    ctx.arc(0, 0, halfsize, 0, 2*Math.PI, true);
    ctx.fill();
    
    ctx.beginPath();
    ctx.arc(0, 0, halfsize, 0, 2*Math.PI, true);
    ctx.clip();
}    

StarMap.prototype.setOrient = function (orient) {
    this.orient = orient;
};

/** A planet object.
 * @constructor
 */
StarMap.Planet = function (pl, size, color) {
    this.pl = pl;
    this.size = size;
    this.color = color;
}

StarMap.Planet.prototype.getCoord = function (jct, earthPos, equ2ecl) {
    var pos = this.pl['keplerCoord'](jct);
    return new StarJs.Vector.Polar3(equ2ecl.apply(pos.sub(earthPos)));
};

/** A Moon object.
 * @constructor
 */
StarMap.Moon = function (size, color) {
    this.size = size;
    this.color = color;
};

StarMap.Moon.prototype.pl = { name: 'Moon' };

StarMap.Moon.prototype.getCoord = function (jct, earthPos, equ2ecl) {
    // earthPos and equ2ecl are ignored
    var pos = StarJs.Solar.approxMoon(jct);
    return {'phi': pos['ra'], 'theta': pos['dec']};
};

StarMap.PLANETS = [
    new StarMap.Planet(StarJs.Solar.BODIES['Sun'], 20, '#FF0'),
    new StarMap.Moon(20, '#880'),
    new StarMap.Planet(StarJs.Solar.BODIES['Mercury'], 3, '#888'),
    new StarMap.Planet(StarJs.Solar.BODIES['Venus'], 4, '#AAA'),
    new StarMap.Planet(StarJs.Solar.BODIES['Mars'], 4, '#F80'),
    new StarMap.Planet(StarJs.Solar.BODIES['Jupiter'], 6, '#FB0'),
    new StarMap.Planet(StarJs.Solar.BODIES['Saturn'], 6, '#AA0'),
    new StarMap.Planet(StarJs.Solar.BODIES['Uranus'], 6, '#CAF'),
    new StarMap.Planet(StarJs.Solar.BODIES['Neptune'], 6, '#CAF')
];

StarMap.EARTH = StarJs.Solar.BODIES['Earth'];

/** Celestial path (of a comet, asteroid, satellite etc.).
 * @constructor
 */
StarMap.Path = function (labels, ras, des, style, labelStyle) {
    this.labels = labels;
    this.ras = ras;
    this.des = des;
    this.style = style;
    this.labelStyle = labelStyle;
};

StarMap.Path.prototype.draw = function (ctx, proj) {
    var labels = this.labels, ra = this.ras, de = this.des;
    var pts = Array(de.length);
    
    for (var key in this.style) {
        ctx[key] = this.style[key];
    }

    ctx.beginPath();
    for (var i = 0; i < de.length; ++i) {
        var cp = pts[i] = proj.projectObj(Math.PI*de[i]/180.0,
                                      Math.PI*ra[i]/12.0);
        if (cp[2]) {
            var xx = cp[0];
            var yy = cp[1];
            ctx.lineTo(xx, yy);
        }
    }
    ctx.stroke();

    for (key in this.labelStyle) {
        ctx[key] = this.labelStyle[key];
    }
    for (i = 0; i < de.length; ++i) {
        var cp = pts[i];
        if (cp[2]) {
            xx = cp[0];
            yy = cp[1];
            ctx.beginPath();
            ctx.arc(xx, yy, 2, 0, 2*Math.PI, true);
            ctx.fill();
            if (ctx.fillText) {
                ctx.fillText(labels[i], xx, yy - 4);
            }
        }
    }

};

/** Telrad pattern.
 * @constructor
 */
StarMap.Telrad = function (lat, lon) {
    this.lat = lat;
    this.lon = lon;
};

StarMap.Telrad.prototype.draw = function (ctx, proj) {
    var p = proj.projectObj(this.lat, this.lon);
    var g05 = proj.projectCircle2(p, 0.5/180*Math.PI);
    var g20 = proj.projectCircle2(p, 2.0/180*Math.PI);
    var g40 = proj.projectCircle2(p, 4.0/180*Math.PI);

    var h = Math.floor(this.size/2);

    ctx.strokeStyle = 'rgba(255,0,0,0.6)';
    ctx.lineWidth = g40.rad/9;

    function drawBullEye(g) {
        var D = 0.15;
        ctx.beginPath();
        ctx.arc(g.x, g.y, g.rad, D, 0.5*Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(g.x, g.y, g.rad, 0.5*Math.PI+D, Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(g.x, g.y, g.rad, Math.PI+D, 1.5*Math.PI-D, false);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(g.x, g.y, g.rad, 1.5*Math.PI+D, 2*Math.PI-D, false);
        ctx.stroke();
    }

    ctx.beginPath();
    ctx.arc(g05.x, g05.y, g05.rad, 0, 2*Math.PI, true);
    ctx.stroke();

    drawBullEye(g20);
    drawBullEye(g40);
};

/** Constellation boundaries.
 * @constructor
 */
StarMap.ConstellationBoundaries = function (boundaries, epoch) {
    var DEG2RAD = StarJs.Math.DEG2RAD;

    if (typeof epoch === 'undefined') {
        epoch = 0; // J2000
    }
    var prec = StarJs.Coord.precessionEquMatrix((1875-2000)/100.0, epoch);

    // Precess polar point
    var polar = new StarJs.Vector.Vector3(0,0,1);
    this.polarPrec = new StarJs.Vector.Polar3(prec.apply(polar));

    // Precess point
    var len = boundaries.length;
    var result = Array(len);
    for (var i = 0; i < len; ++i) {
        var pt = boundaries[i].slice(0);
        var v = new StarJs.Vector.Polar3(15*Math.PI*pt[0]/180,
                                         DEG2RAD*pt[1]).toVector3();
        var p = new StarJs.Vector.Polar3(prec.apply(v));
        pt[0] *= 15;
        pt.push(p['phi']/DEG2RAD)
        pt.push(p['theta']/DEG2RAD);
        result[i] = pt;
    }

    this.boundaries = result;
};

StarMap.ConstellationBoundaries.prototype.draw = function (ctx, proj) {
    var DEG2RAD = StarJs.Math.DEG2RAD;
    function angSep(a1, a2) {
        return StarJs.Math.mod(a1-a2, 2 * Math.PI);
    };
    // Constellation boundaries
    var cstn = null;
    // TODO: make a = 0.5 after boundaries line merging.  Now many
    // lines are drawn twice, and they color differs (somewhat
    // brighter) if alpha is used.
    ctx.strokeStyle = 'rgba(128,0,128,1)';
    var prev = null, first, polarPrec = this.polarPrec;
    function drawCnstBnd(prev, l) {
	var seg, flip, a1, a2;
	if (prev[1] === l[1]) {
            a1 = DEG2RAD*prev[0];
            a2 = DEG2RAD*l[0];
	    seg = proj.projectSegment(
                // Center
                polarPrec['theta'], polarPrec['phi'],
                // Radius
                Math.PI/2-DEG2RAD*l[1],
                // Point 1
                prev[4]*DEG2RAD,
                prev[3]*DEG2RAD,
                // Point 2
                l[4]*DEG2RAD,
                l[3]*DEG2RAD
	    );
	    flip = angSep(a1, a2) < angSep(a2, a1) !== seg.flip;
	    gr = false;
	} else {
            a1 = DEG2RAD*prev[1];
            a2 = DEG2RAD*l[1];
            seg = proj.projectGreatSegment(
                prev[4]*DEG2RAD,
                prev[3]*DEG2RAD,
                // Point 2
                l[4]*DEG2RAD,
                l[3]*DEG2RAD
            );
	    gr = true;
	    flip = seg.flip;
	}

	ctx.beginPath();
	switch (seg.type) {
	case 'circle':
	    ctx.arc(seg.x, seg.y,
		    seg.r,
		    seg.a1, seg.a2,
		    flip);
	    break;
	case 'line':
	    // TODO sometimes lines shouldn't be drawn if their
	    // central point pass through infinity, or at least
	    // they should be drawn more intelligently.  Can be
	    // such lines visible in viewport?
	    ctx.moveTo(seg.x1, seg.y1);
	    ctx.lineTo(seg.x2, seg.y2);
	    break;
	}
	ctx.stroke();
    }
    var boundaries = this.boundaries;
    for (var j = 0; j < boundaries.length; ++j) {
        var l = boundaries[j], a1, a2, gr;
        if (cstn === l[2]) {
	    drawCnstBnd(prev, l);
        } else {
            if (prev !== null) {
                drawCnstBnd(prev, first);
            }
            first = prev = l;
            cstn = l[2];
        }
        prev = l;
    }
};

/** An object catalogue (Messier, Caldwell, NGC, etc).
 * @constructor
 */
StarMap.Catalogue = function (name, data, colors, renderer) {
    function messierColor(mag) {
        var v = Math.min(15, Math.floor(19-mag));
        var h = v.toString(16);
        return '#'+h+h+h;
    }

    this.name = name;
    this.data = data.slice(0);
    this.colors = colors;
    
    this.renderer = renderer || function (ctx, cm, cc, colors) {
        ctx.beginPath();
        ctx.strokeStyle = colors[cc[2]] || messierColor(cc[5]);
        ctx.arc(cm[0], cm[1], 4, 0, 2*Math.PI, true);
        ctx.stroke();
    };
};

StarMap.Catalogue.prototype.draw = function (ctx, proj) {
    var data = this.data;
    var len = data.length, cc, cm;

    for (var i = 0; i < len; ++i) {
        cc = data[i];
        cm = proj.projectObj(cc[4], 15*cc[3]);
        if (cm[2]) {
            this.renderer(ctx, cm, cc, this.colors);
        }
    }
};

/** A graticule.
 * @constructor
 */
StarMap.Graticule = function (lon) {
    this.lon = lon;
};

StarMap.Graticule.prototype.draw = function (ctx, proj) {
    ctx.strokeStyle = '#448';
    var halfsize = Math.floor(this.size / 2);
    for (var i = -80; i < 90; i += 10) {
        var p = proj.projectParallel(Math.PI*i/180);
        ctx.beginPath();
        ctx.lineWidth = (i === 0) ? 1.7 : 1;
        switch (p.type) {
        case 'line':
            ctx.moveTo(p.x-halfsize*p.vx,
                       p.y-halfsize*p.vy);
            ctx.lineTo(p.x+halfsize*p.vx,
                       p.y+halfsize*p.vy);
            break;
        case 'circle':
            ctx.arc(p.x, p.y, p.r, 0, 2*Math.PI, true);
            break;
        }
        ctx.stroke();
    }
    for (i = 0; i < 180; i += 15) {
        var p = proj.projectMeridian(Math.PI*i/180);
        ctx.beginPath();
        ctx.lineWidth = (i === 0 || i === -180) ? 1.7 : 1;
        switch (p.type) {
        case 'line':
            ctx.moveTo(p.x-halfsize*p.vx,
                       p.y-halfsize*p.vy);
            ctx.lineTo(p.x+halfsize*p.vx,
                       p.y+halfsize*p.vy);
            break;
        case 'circle':
            ctx.arc(p.x, p.y, p.r,
                    0, 2*Math.PI, true);
            break;
        }
        ctx.stroke();
    }

    ctx.lineWidth = 2;
    ctx.beginPath();
    p = proj.projectParallel(Math.PI/2-this.lon);
    switch (p.type) {
    case 'line':
        ctx.moveTo(p.x-halfsize*p.vx,
                   p.y-halfsize*p.vy);
        ctx.lineTo(p.x+halfsize*p.vx,
                   p.y+halfsize*p.vy);
        break;
    case 'circle':
        ctx.arc(p.x, p.y, p.r, 0, 2*Math.PI, true);
        break;
    }
    ctx.stroke();
    ctx.beginPath();
    p = proj.projectParallel(-Math.PI/2-this.lon);
    switch (p.type) {
    case 'line':
        ctx.moveTo(p.x-halfsize*p.vx,
                   p.y-halfsize*p.vy);
        ctx.lineTo(p.x+halfsize*p.vx,
                   p.y+halfsize*p.vy);
        break;
    case 'circle':
        ctx.arc(p.x, p.y, p.r, 0, 2*Math.PI, true);
        break;
    }
    ctx.stroke();
    ctx.lineWidth = 1;
};

/** Ecliptics circle.
 * @constructor
 */
StarMap.Ecliptics = function () {
};

StarMap.Ecliptics.prototype.draw = function (ctx, proj) {
    var eclp = proj.projectCircle(Math.PI/2 + StarJs.Solar.EPS, Math.PI/2, Math.PI/2);
    if (eclp.type === 'circle') {
        ctx.beginPath();
        ctx.strokeStyle = 'yellow';
        ctx.arc(eclp.x, eclp.y, eclp.rad, 0, 2*Math.PI, true);
        ctx.stroke();
    }
};

/** Any object.
 * @constructor
 */
StarMap.Object = function (params, color, mjd, jct) {
    this.params = params;
    this.color = color;
    this.pqr = StarJs.Kepler.gaussVec(params.node, params.incl, params.peri);
    this.mjd = mjd;
    this.jct = jct;
};

StarMap.Object.prototype.draw = function (ctx, proj) {
    var pos = StarJs.Kepler.keplerPos(StarJs.Solar.GM,
                                      this.params.t0,
                                      this.mjd,
                                      this.params.q,
                                      this.params.e,
                                      this.pqr);
    var earthPos = StarMap.EARTH['keplerCoord'](this.jct);
    var ecl2equ = StarJs.Coord.ecl2equMatrix(this.jct);
    pos = new StarJs.Vector.Polar3(ecl2equ.apply(pos.sub(earthPos)));
    var cm = proj.projectObj(pos['theta'], pos['phi']);

    ctx.strokeStyle = ctx.fillStyle = this.color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(cm[0], cm[1], 6, 0, 2*Math.PI, true);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(cm[0], cm[1], 4, 0, 2*Math.PI, true);
    ctx.fill();

    if (this.params.title) {
        ctx.fillText(this.params.title, cm[0]+8, cm[1]+8);
    }
};

StarMap.prototype.setPos = function (lat, lon, time) {
    if (typeof time === 'undefined') {
        time = +new Date();
    } else if (typeof time !== 'number') {
        time = +time;
    }

    this.lat = lat;
    this.lon = lon;
    this.time = time;
};

StarMap.prototype.draw = function () {
    var lat = this.lat, lon = this.lon, time = this.time;

    var Ti = StarJs.Time;

    var mjd = Ti.time2mjd(time);
    var gms_t = Ti.gmst(mjd);

    /** @const */
    var DEG2RAD = StarJs.Math.DEG2RAD;
    lat *= DEG2RAD;
    lon *= DEG2RAD;

    lat += gms_t;

    this['proj'].setCoords(lon, lat);

    var ortho = this['proj'].projectPoints(this.stars);
    var cst = [], i, j, slen = ortho.length, co = this.cnstltns, clen = co.length, halfsize = Math.floor(this.size/2);
    
    var ctx = this.ctx;
    ctx.save();

    ctx.translate(halfsize, halfsize);
    ctx.rotate(this.orient);

    this.drawBg();

    // Draw graticule
    (new StarMap.Graticule(lon)).draw(ctx, this['proj']);


    // Boundaries
    (new StarMap.ConstellationBoundaries(window.CON_BOUND_18)).draw(ctx, this['proj']);

    // Constellations
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(255,255,255,0.7)';
    for (j = clen; j--; ) {
        var s = co[j][0], e = co[j][1];
        var so = ortho[s], eo = ortho[e];
        if (so[3] || eo[3]) {
            ctx.moveTo(so[1], so[2]);
            ctx.lineTo(eo[1], eo[2]);
        }
    }
    ctx.stroke();

    // Draw ecliptics
    (new StarMap.Ecliptics()).draw(ctx, this['proj']);

    // Stars
    ctx.fillStyle = '#FFF';
    for (i = 0; i < slen; ++i) {
        var s = ortho[i];
        if (s[3]) {
            ctx.beginPath();
            ctx.arc(s[1], s[2],
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
    if (this.prop && this.prop['messier']) {
        (new StarMap.Catalogue("Messier", this.prop['messier'], this.prop['messier_colors'])).draw(ctx, this['proj']);
    }
    // Draw Caldwell objects
    if (this.prop && this.prop['caldwell']) {
        (new StarMap.Catalogue("Caldwell", this.prop['caldwell'], this.prop['caldwell_colors'])).draw(ctx, this['proj']);
    }

    // Draw planets
    var jct = Ti.mjd2jct(mjd);
    if (this.planets) {
        var earthPos = StarMap.EARTH['keplerCoord'](jct);
        var equ2ecl = StarJs.Coord.ecl2equMatrix(jct);
        for (i = 0; i < StarMap.PLANETS.length; ++i) {
            var planet = StarMap.PLANETS[i];
            var cc = planet.getCoord(jct, earthPos, equ2ecl);
            
            var cm = this['proj'].projectObj(cc['theta'], cc['phi']);
            if (cm[2]) {
                ctx.beginPath();
                ctx.fillStyle = planet.color;
                var xx = cm[0], yy = cm[1];
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
    (new StarMap.Telrad(0.901, 0.451)).draw(ctx, this['proj']);
    (new StarMap.Telrad(0.301, 2.151)).draw(ctx, this['proj']);

    /*
Ephemerides from http://www.minorplanetcenter.net/iau/MPEph/MPEph.html

C/2009 P1 (Garradd)

Perturbed ephemeris below is based on elements from MPEC 2012-A52.

   CK09P010
Date       UT      R.A. (J2000) Decl.    Delta     r     El.    Ph.   m1     Sky Motion
            h m s                                                            "/min    P.A.
2012 01 17 000000   17.4533  +32.839     1.749   1.586   64.0  33.9   7.2    1.17    350.9
2012 01 22 000000   17.4175  +35.301     1.683   1.602   68.1  34.8   7.2    1.34    348.6
2012 01 27 000000   17.3649  +38.107     1.615   1.621   72.6  35.4   7.1    1.54    346.0
2012 02 01 000000   17.2890  +41.290     1.548   1.642   77.3  35.8   7.1    1.77    343.2
2012 02 06 000000   17.1802  +44.877     1.482   1.666   82.2  35.9   7.1    2.02    340.1
2012 02 11 000000   17.0241  +48.874     1.421   1.692   87.3  35.6   7.0    2.29    336.3
2012 02 16 000000   16.7981  +53.248     1.367   1.720   92.4  35.0   7.0    2.57    331.7
2012 02 21 000000   16.4653  +57.889     1.321   1.750   97.4  34.1   7.0    2.85    325.6
2012 02 26 000000   15.9670  +62.547     1.288   1.782  102.1  32.9   7.1    3.09    317.2
2012 03 02 000000   15.2189  +66.752     1.269   1.815  106.2  31.6   7.1    3.26    305.0
2012 03 07 000000   14.1497  +69.747     1.267   1.851  109.4  30.4   7.2    3.35    288.1
2012 03 12 000000   12.8391  +70.676     1.283   1.888  111.4  29.4   7.3    3.33    267.6
2012 03 17 000000   11.5897  +69.266     1.317   1.926  112.0  28.6   7.4    3.20    247.8
2012 03 22 000000   10.6405  +66.145     1.368   1.965  111.4  28.2   7.6    3.00    232.4
2012 03 27 000000   09.9973  +62.191     1.436   2.005  109.7  27.9   7.8    2.75    221.2
2012 04 01 000000   09.5739  +58.016     1.517   2.047  107.0  27.8   8.0    2.48    213.1
2012 04 06 000000   09.2944  +53.951     1.611   2.089  103.8  27.7   8.2    2.21    206.8
2012 04 11 000000   09.1087  +50.147     1.714   2.132  100.1  27.6   8.5    1.96    201.6
2012 04 16 000000   08.9858  +46.657     1.825   2.176   96.2  27.3   8.7    1.74    196.9
2012 04 21 000000   08.9067  +43.487     1.943   2.221   92.2  26.9   8.9    1.54    192.6
2012 04 26 000000   08.8588  +40.618     2.065   2.266   88.0  26.4   9.1    1.38    188.4
2012 05 01 000000   08.8340  +38.022     2.190   2.311   83.9  25.7   9.3    1.24    184.2
2012 05 06 000000   08.8263  +35.669     2.317   2.358   79.8  24.9   9.5    1.12    180.2
2012 05 11 000000   08.8317  +33.526     2.446   2.404   75.7  24.0   9.8    1.03    176.2
2012 05 16 000000   08.8474  +31.567     2.574   2.451   71.6  23.0   9.9    0.95    172.2
2012 05 21 000000   08.8713  +29.766     2.703   2.498   67.6  22.0  10.1    0.88    168
    */

    // Draw path of C/2009 R1 (McNaught)
    var C2009P1_DA = ["2012 01 17", "2012 01 22", "2012 01 27",
                      "2012 02 01", "2012 02 06", "2012 02 11",
                      "2012 02 16", "2012 02 21", "2012 02 26",
                      "2012 03 02", "2012 03 07", "2012 03 12",
                      "2012 03 17", "2012 03 22", "2012 03 27",
                      "2012 04 01", "2012 04 06", "2012 04 11"];
    var C2009P1_RA = [
        17.4533, 17.4175, 17.3649,
        17.2890, 17.1802, 17.0241,
        16.7981, 16.4653, 15.9670,
        15.2189, 14.1497, 12.8391,
        11.5897, 10.6405,  9.9973,
         9.5739,  9.2944,  9.1087
    ];
    var C2009P1_DE = [
            +32.839, +35.301, +38.107,
            +41.290, +44.877, +48.874,
            +53.248, +57.889, +62.547,
            +66.752, +69.747, +70.676,
            +69.266, +66.145, +62.191,
            +58.016, +53.951, +50.147
    ];
    var C2009P1 = new StarMap.Path(C2009P1_DA, C2009P1_RA, C2009P1_DE, {
        'fillStyle': 'green',
        'strokeStyle': 'green',
        'lineWidth': 1
    }, {
        'fillStyle': '#8F8',
        'strokeStyle': '#8F8'
    });
    C2009P1.draw(ctx, this['proj']);
    
    /*
    http://www.minorplanetcenter.net/mpec/K12/K12A52.html

    C/2009 P1 (Garradd)
Epoch 2011 Dec. 25.0 TT = JDT 2455920.5
T 2011 Dec. 23.67758 TT                                 MPC
q   1.5505367            (2000.0)            P               Q
z  -0.0006782      Peri.   90.74770     -0.16661319     -0.82691127
 +/-0.0000002      Node   325.99770     -0.58719595     +0.52078702
e   1.0010516      Incl.  106.17749     +0.79211171     +0.21212879
From 4943 observations 2009 Aug. 13-2011 Dec. 4, mean residual 0".4.
    */
    // Draw current position of the comet
    var C2009P1_param = {
        title: 'C/2009 P1 (Garradd)',
        t0: 2455920.5-2400000.5,
        q: 1.5505367,
        z: -0.0006782,
        e: 1.0010516,
        peri: 90.74770*DEG2RAD,
        node: 325.99770*DEG2RAD,
        incl: 106.17749*DEG2RAD
    };

    (new StarMap.Object(C2009P1_param, 'rgba(220, 255, 220, 0.7)', mjd, jct)).draw(ctx, this['proj']);


    var lutetia = {
        title: 'Lutetia',
        t0: 2455400.5-2400000.5,
        q: 2.039175887090527,
        e: 0.1628669085598194,
        peri: 250.192513362607*DEG2RAD,
        node: 80.89961160386014*DEG2RAD,
        incl: 3.063753824680438*DEG2RAD
    };

    (new StarMap.Object(lutetia, 'rgba(255, 255, 255, 0.7)', mjd, jct)).draw(ctx, this['proj']);

    
    // Draw sides of Earth
    if (ctx.fillText) {
        //ctx.translate(halfsize, halfsize);
    
        ctx.fillStyle = 'gold';
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 2;
        ctx.font = '20px serif';

        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.strokeText("S", 0, halfsize-5);
        ctx.fillText("S", 0, halfsize-5);

        ctx.rotate(Math.PI/2+0.00001/* Opera 10.60 suxx */);
        ctx.strokeText("E", 0, halfsize-5);
        ctx.fillText("E", 0, halfsize-5);

        ctx.rotate(Math.PI/2);
        ctx.strokeText("N", 0, halfsize-5);
        ctx.fillText("N", 0, halfsize-5);

        ctx.rotate(Math.PI/2);
        ctx.strokeText("W", 0, halfsize-5);
        ctx.fillText("W", 0, halfsize-5);
    }
    ctx.restore();
};

window['StarMap']=StarMap;
StarMap.prototype['setPos'] = StarMap.prototype.setPos;
StarMap.prototype['draw'] = StarMap.prototype.draw;
