var StarJs = {};
StarJs.Math = {};

/**
 * Scale factor for converting degrees to radians.
 * @const
 * @type {number}
 */
StarJs.Math.DEG2RAD = Math.PI / 180.0;
/**
 * Scale factor for converting radians to degrees.
 * @const
 * @type {number}
 */
StarJs.Math.RAD2DEG = 180.0 / Math.PI;
/**
 * Scale factor for converting radians to arcseconds.
 * @const
 * @type {number}
 */
StarJs.Math.RAD2ARCS = 648000.0 / Math.PI;
/**
 * Scale factor for converting arcseconds to radians.
 * @const
 * @type {number}
 */
StarJs.Math.ARCS2RAD = Math.PI / 648000.0;
/**
 * 2*PI
 * @const
 * @type {number}
 */
StarJs.Math.PI2 = 2.0 * Math.PI;

/**
 * Scale factor for converting radians to hour measure.
 * @const
 * @type {number}
 */
StarJs.Math.ANGLE2HMS = 12.0 / Math.PI;

/**
 * Return square of x (x*x)
 * @param {number} x argument
 * @return {number} x*x
 */
StarJs.Math.sqr = function (x) {
    return x * x;
};

/**
 * Return fractional part of argument
 * @param {number} x argument
 * @return {number} fractioal part
 */
StarJs.Math.frac = function (x) {
    return x - Math.floor(x);
};

/**
 * Modulo
 * @param {number} x dividend
 * @param {number} r divisor
 * @return {number} remainder
 */
StarJs.Math.mod = function (x, r) {
    return r * StarJs.Math.frac(x / r);
};

/**
 * Degree-minute-second object.
 * Either 1 or 4 arguments are accepted.
 * @constructor
 * @param {number} sign_or_angle sign if four values are passed, angle
 *                 in degrees otherwise
 * @param {number=} degree Integer degree part (optional)
 * @param {number=} minute Integer minute part (optional)
 * @param {number=} second Seconds (optional)
 */
StarJs.Math.Dms = function (sign_or_angle, degree, minute, second) {
    if (arguments.length === 4) {
        this['sign'] = sign_or_angle;
        this['degree'] = degree;
        this['minute'] = minute;
        this['second'] = second;
    } else {
        var angle = sign_or_angle;
        this['sign'] = 1;
        if (angle < 0) {
            angle = -angle;
            this['sign'] = -1;
        }
        this['degree'] = Math.floor(angle);
        angle = 60 * (angle - this['degree']);
        this['minute'] = Math.floor(angle);
        angle = 60 * (angle - this['minute']);
        this['second'] = angle;
    }
};

/**
 * Convert angle DMS (degree-minute-second) to float degree value.
 * @param {number} sign
 * @param {number} deg
 * @param {number} minute
 * @param {number} second
 */
StarJs.Math.dms2deg = function (sign, deg, minute, second) {
    return sign * (deg + minute / 60.0 + second / 3600.0);
};

/**
 * Convert angle (degree-minute-second) to float degree value.
 */
StarJs.Math.Dms.prototype.dms2deg = function () {
    return StarJs.Math.dms2deg(this['sign'], this['degree'], this['minute'], this['second']);
};

/**
 * Convert float degree value to DMS (degree-minute-second).
 * @param {number} deg Angle
 */
StarJs.Math.deg2dms = function (deg) {
    return new StarJs.Math.Dms(deg);
};

/**
 * Convert radians to hour measure.
 * @param {number} angle Angle in radians.
 */
StarJs.Math.angle2hms = function (angle) {
    angle = StarJs.Math.mod(angle, StarJs.Math.PI2);
    var a = StarJs.Math.ANGLE2HMS * angle, res = StarJs.Math.deg2dms(a);
    // Change field names and remove sign field as it is always 1.
    res['hour'] = res['deg'];
    delete res['deg'];
    delete res['sign'];
    
    return res;
};

/**
 * @param {number} ym
 * @param {number} y0
 * @param {number} yp
 */
StarJs.Math.quadInterpolation = function (ym, y0, yp) {
    var a = 0.5 * (yp + ym) - y0, b = 0.5 * (yp - ym), c = y0, xe = -b / (2 * a), ye = (a * xe + b) * xe + c;
    var dis = b * b - 4 * a * c, roots = [], dx, r1, r2;
    if (dis >= 0) {
        dx = 0.5 * Math.sqrt(dis) / Math.abs(a);

        r1 = xe - dx;
        r2 = xe + dx;
        
        if (Math.abs(r1) <= 1.0) {
            roots.push(r1);
        }
        if (Math.abs(r2) <= 1.0) {
            roots.push(r2);
        }
    }
    return {'xe': xe, 'ye': ye, 'roots': roots};
};

/**
 * Hyperbolic sinus.
 * @param {number} x
 * @return {number}
 */
StarJs.Math.sinh = function (x) {
    return (Math.exp(x) - Math.exp(-x))/2;
};

/**
 * Hyperbolic cosine.
 * @param {number} x
 * @return {number}
 */
StarJs.Math.cosh = function (x) {
    return (Math.exp(x) + Math.exp(-x))/2;
};

StarJs.Vector = {};

/** @constructor
 * @param {number} x
 * @param {number} y
 * @param {number} z
 */
StarJs.Vector.Vector3 = function (x, y, z) {
    this['x'] = x;
    this['y'] = y;
    this['z'] = z;
};

(function (p) {
    p.len = function () {
        return Math.sqrt(this['x'] * this['x'] + this['y'] * this['y'] + this['z'] * this['z']);
    };

    /** Return new Vector3 as sum of this and the v3.
     * @param {StarJs.Vector.Vector3} v3
     */
    p.add = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] + v3['x'],
                                         this['y'] + v3['y'],
                                         this['z'] + v3['z']);
    };

    /** Return new Vector3 as v3 substracted from this.
     * @param {StarJs.Vector.Vector3} v3
     */
    p.sub = function (v3) {
        return new StarJs.Vector.Vector3(this['x'] - v3['x'],
                                         this['y'] - v3['y'],
                                         this['z'] - v3['z']);
    };

    p.neg = function () {
        return new StarJs.Vector.Vector3(-this['x'], -this['y'], -this['z']);
    };

    /** Return new Vector3 multiplied by number a.
     * @param {number} a
     */
    p.scale = function (a) {
        return new StarJs.Vector.Vector3(a * this['x'], a * this['y'], a * this['z']);
    };

    /** Copy this vector.
     */
    p.clone = function () {
        return new StarJs.Vector.Vector3(this['x'], this['y'], this['z']);
    };
}(StarJs.Vector.Vector3.prototype));

/** @constructor
 * @param {number|StarJs.Vector.Vector3} az_or_v3
 * @param {number=} elev
 * @param {number=} rad
 */
StarJs.Vector.Polar3 = function (az_or_v3, elev, rad) {
    var alen = arguments.length;
    if (alen === 2) {
        rad = 1.0;
    }
    if (alen === 2 || alen === 3) {
        this['phi'] = az_or_v3;
        this['theta'] = elev;
        this['rad'] = rad;
    } else {
        var v3 = az_or_v3;
        var rho2 = v3['x'] * v3['x'] + v3['y'] * v3['y'];
        this['rad'] = Math.sqrt(rho2 + v3['z'] * v3['z']);
        this['phi'] = (v3['x'] === 0.0 && v3['y'] === 0.0) ?  0.0 : Math.atan2(v3['y'], v3['x']);
        if (this['phi'] < 0.0) {
            this['phi'] += StarJs.Math.PI2;
        }
        var rho = Math.sqrt(rho2);
        this['theta'] = (v3['z'] === 0.0 && rho === 0.0) ? 0.0 : Math.atan2(v3['z'], rho);
    }
};

(function (p) {
    p.toVector3 = function () {
        var ct = Math.cos(this['theta']);
        var rad = this['rad']
        var x = rad * ct * Math.cos(this['phi']);
        var y = rad * ct * Math.sin(this['phi']);
        var z = rad * Math.sin(this['theta']);
        return new StarJs.Vector.Vector3(x, y, z);
    };

    p.clone = function () {
        return new StarJs.Vector.Polar3(this['rad'], this['phi'], this['theta']);
    };
}(StarJs.Vector.Polar3.prototype));

/** @constructor
 * @param {(StarJs.Vector.Vector3|boolean)=} v2_or_reuse_flag
 * @param {StarJs.Vector.Vector3=} v3
 */
StarJs.Vector.Matrix3 = function (v1, v2_or_reuse_flag, v3) {
    if (arguments.length === 3) {
        var v2 = v2_or_reuse_flag;  // Just for brevity.
        this['mat'] = [[v1['x'], v1['y'], v1['z']],
                       [v2['x'], v2['y'], v2['z']],
                       [v3['x'], v3['y'], v3['z']]];
    } else {
        if (v2_or_reuse_flag) {
            this['mat'] = v1;
        } else {
            this['mat'] = [[v1[0][0], v1[0][1], v1[0][2]],
                        [v1[1][0], v1[1][1], v1[1][2]],
                        [v1[2][0], v1[2][1], v1[2][2]]];
        }
    }
};

(function (p) {
    /** Multiply matrix-to-vector multiplication.
     * @param {StarJs.Vector.Vector3} v
     */
    p.apply = function (v) {
        var l = this['mat'][0];
        var x = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        l = this['mat'][1];
        var y = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        l = this['mat'][2];
        var z = l[0] * v['x'] + l[1] * v['y'] + l[2] * v['z'];

        return new StarJs.Vector.Vector3(x, y, z);
    };

    /** Perform a matrix multiplication, returning a new matrix.
     * @param {StarJs.Vector.Matrix3} matrix right matrix.
     */
    p.mult = function (matrix) {
        var res = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        var mat = matrix['mat'];
        for (var i = 0; i < 3; i += 1) {
            var tline = this['mat'][i];
            var rline = res[i];
            for (var j = 0; j < 3; j += 1) {
                for (var k = 0; k < 3; k += 1) {
                    rline[j] += tline[k] * mat[k][j];
                }
            }
        }
        return new StarJs.Vector.Matrix3(res, true);
    };
}(StarJs.Vector.Matrix3.prototype));

/** Return x-rotation matrix by angle phi.
 * @param {number} phi.
 */
StarJs.Vector.Matrix3.r_x = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[1.0, 0.0, 0.0],
                                      [0.0,  cp,  sp],
                                      [0.0, -sp,  cp]]);
};

/** Return y-rotation matrix by angle phi.
 * @param {number} phi.
 */
StarJs.Vector.Matrix3.r_y = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp, 0.0, -sp],
                                      [0.0, 1.0, 0.0],
                                      [ sp, 0.0,  cp]]);
};

/** Return z-rotation matrix by angle phi.
 * @param {number} phi.
 */
StarJs.Vector.Matrix3.r_z = function (phi) {
    var cp = Math.cos(phi), sp = Math.sin(phi);
    return new StarJs.Vector.Matrix3([[ cp,  sp, 0.0],
                                      [-sp,  cp, 0.0],
                                      [0.0, 0.0, 1.0]]);
    
};
StarJs.Time = {};

StarJs.Time['DEFAULT_JULIAN_DATE'] = {'year': 1582, 'month': 10, 'day': 4};
StarJs.Time['DEFAULT_JULIAN_JD'] = 2299161;
StarJs.Time['JD_MJD'] = 2400000.5;

/** Convert Date or timestamp to mjd float.
 * @param {Date|number} t time
 */
StarJs.Time.time2mjd = function (t) {
    if (typeof t !== 'number') {
        t = t.getTime();
    }
    return t / 86400000.0 + 40587.0;
};

/** Convert mjd to date and time object, taking into account Julian
 * calendar.
 * @param jul {number=} JD of Julian date change.
 */
StarJs.Time.mjd2dt = function (mjd, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time['DEFAULT_JULIAN_JD'];
    }
    var a = Math.floor(mjd) + 2400001, b, c, d, e, f, day, mon, year, hour, t;
    if (a < jul) {
        // Julian calendar
        b = 0;
        c = a + 1524;
    } else {
        // Grigorian calendar
        b = Math.floor((a - 1867216.25) / 36524.25);
        c = a + b - Math.floor(b / 4) + 1525;
    }
    d = Math.floor((c - 122.1) / 365.25);
    e = 365 * d + Math.floor(d / 4);
    f = Math.floor((c - e) / 30.6001);
    day = c - e - Math.floor(30.6001 * f);
    mon = f - 1 - 12 * Math.floor(f / 14);
    year = d - 4715 - Math.floor((7 + mon) / 10);
    hour = 24.0 * (mjd - Math.floor(mjd));
    t = StarJs.Time.hour2hms(hour);
    t['year'] = year;
    t['month'] = mon;
    t['day'] = day;
    return t;
};

/** (h, m, s) to floating-point hour.
 * @param {number} h hour
 * @param {number} m minute
 * @param {number} s second
 */
StarJs.Time.hms2hour = function (h, m, s) {
    return h + (m / 60.0) + (s / 3600.0);
};

/** Convert floating-point hour to hours, minutes and seconds,
 * returing an object.
 * @param {number} hour
 */
StarJs.Time.hour2hms = function (hour) {
    var dms = new StarJs.Math.Dms(hour);
    return {'hour': dms['degree'], 'minute': dms['minute'], 'second': dms['second']};
};

/** Convert data-time object to MDJ float, taking into account a Julian
 * calendar.
 * @param {number} dt
 * @param {number=} jul  JD of Julian date change.
 */
StarJs.Time.dt2mjd = function (dt, jul) {
    if (typeof jul === 'undefined') {
        jul = StarJs.Time['DEFAULT_JULIAN_DATE'];
    }

    var year = dt['year'], mon = dt['month'], b;
    if (mon <= 2) {
        mon += 12;
        year -= 1;
    }
    if (year <= jul['year'] && mon <= jul['month'] && dt['day'] <= jul['day']) {
        // Julian
        b = -2 + Math.floor((year + 4716) / 4) - 1179;
    } else {
        // Gregorian
        b = Math.floor(year / 400) - Math.floor(year / 100) + Math.floor(year / 4);
    }
    var mjdMidnight = 365 * year - 679004 + b + Math.floor(30.6001 * (mon + 1)) + dt['day'];
    var frac = StarJs.Time.hms2hour(dt['hour'], dt['minute'], dt['second']) / 24.0;
    return mjdMidnight + frac;
};

/** Convert MJD to JCT.
 * @param {number} mjd
 */
StarJs.Time.mjd2jct = function (mjd) {
    return (mjd - 51544.5) / 36525.0;
};

/** Convert MJD to Greenwich Mean Sidereal Time float.
 * @param mjd {number}
 */
StarJs.Time.gmst = function (mjd) {
    /* TODO: move to global */
    var SECS = 86400; // 24*60*60 -- number of seconds in day;
    var mjd0 = Math.floor(mjd), ut = SECS * (mjd - mjd0), t0, t, gmst;
    t0 = StarJs.Time.mjd2jct(mjd0);
    t = StarJs.Time.mjd2jct(mjd);
    gmst = 24110.54841 + 8640184.812866 * t0 + 1.0027379093 * ut +
	(0.093104 - (6.2e-6) * t) * t * t;
    return StarJs.Math.PI2 / SECS * StarJs.Math.mod(gmst, SECS);
};
StarJs.Coord = {};

/** Precession matrix (in ecliptical coordinates) from epoch t1 to
 *  epoch t2.
 * @param {number} t1
 * @param {number} t2
 */
StarJs.Coord.precessionEclMatrix = function (t1, t2) {
    var dt = t2 - t1, p1, p2, pa;
    p1 = 174.876383889 * StarJs.Math.DEG2RAD +
	(((3289.4789 + 0.60622 * t1) * t1) +
	 ((-869.8089 - 0.50491 * t1) + 0.03536 * dt) * dt) * StarJs.Math.ARCS2RAD;
    p2 = ((47.0029 - (0.06603 - 0.000598 * t1) * t1) +
	  ((-0.03302 + 0.00598 * t1) + 0.000060 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    pa = ((5029.0966 + (2.22226 - 0.000042 * t1) * t1) +
	  ((1.11113 - 0.000042 * t1) - 0.000006 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-(p1 + pa)).mult(StarJs.Vector.Matrix3.r_x(p2)).mult(StarJs.Vector.Matrix3.r_z(p1));
};

/** Precession matrix (in equatorial coordinates) from epoch t1 to
 *  epoch t2.
 * @param {number} t1
 * @param {number} t2
 */
StarJs.Coord.precessionEquMatrix = function (t1, t2) {
    var dt = t2 - t1, zeta, z, theta;
    zeta = ((2306.2181 + (1.39656 - 0.000139 * t1) * t1) +
	    ((0.30188 - 0.000344 * t1) + 0.017998 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    z = zeta + ((0.79280 + 0.000411 * t1) + 0.000205 * dt) * dt * dt * StarJs.Math.ARCS2RAD;
    theta  = ((2004.3109 - (0.85330 + 0.00217 * t1) * t1) -
	      ((0.42665 + 0.000217 * t1) + 0.041833 * dt) * dt) * dt * StarJs.Math.ARCS2RAD;
    return StarJs.Vector.Matrix3.r_z(-z).mult(StarJs.Vector.Matrix3.r_y(theta)).mult(StarJs.Vector.Matrix3.r_z(-zeta));
    
};

/** Oliquity of the ecliptic.
 *  You may use StarJs.Solar.EPS as constant approximation.
 * @param {number} jct
 */
StarJs.Coord.eclipticObliquity = function (jct) {
    return (23.43929111 - (46.8150 + (0.00059 - 0.001813 * jct) * jct) * jct / 3600.0) * StarJs.Math.DEG2RAD;
};

/** Matrix for conversion from equatorial to ecliptic coordinate system.
 * @param {number} jct
 */
StarJs.Coord.equ2eclMatrix = function (jct) {
    var eps = StarJs.Coord.eclipticObliquity(jct);
    return StarJs.Vector.Matrix3.r_x(eps);
};

/** Matrix for conversion from ecliptic to equatorial coordinate system.
 * @param {number} jct
 */
StarJs.Coord.ecl2equMatrix = function (jct) {
    var eps = StarJs.Coord.eclipticObliquity(jct);
    return StarJs.Vector.Matrix3.r_x(-eps);
};

/** Convert from equatorial to horizontal coordinate system.
 * @param {number} dec
 * @param {number} tau
 * @param {number} lat
 */
StarJs.Coord.equ2hor = function (dec, tau, lat) {
    var hor_vec = StarJs.Vector.Matrix3.r_y(Math.PI / 2 - lat).apply(new StarJs.Vector.Polar3(tau, dec).toVector3());
    var hor_pol = new StarJs.Vector.Polar3(hor_vec);
    return {'h': hor_pol['theta'], 'az': hor_pol['phi']};
};

/** Convert from horizontal to equatorial coordinate system.
 * @param {number} h
 * @param {number} az
 * @param {number} lat
 */
StarJs.Coord.hor2equ = function (h, az, lat) {
    var equ_vec = StarJs.Vector.Matrix3.r_y(-Math.PI / 2 + lat).apply(new StarJs.Vector.Polar3(az, lat).toVector3());
    var equ_pol = new StarJs.Vector.Polar3(equ_vec);
    return {'dec': equ_pol['theta'], 'tau': equ_pol['phi']};
};
StarJs.Kepler = {};

/** Default number of iterations for iterative method. */
StarJs.Kepler['DEFAULT_ITERATIONS'] = 100;
/** Default iteration precision for iterative method. */
StarJs.Kepler['DEFAULT_PRECISION'] = 1e-9;

/** Compute eccentric anomaly for elliptical orbit.
 * @param {number} m Mean anomaly.
 * @param {number} ec Eccentricity (ec < 1)
 * @param {number=} maxiter Optional number of iterations.
 */
StarJs.Kepler.eccAnomaly = function (m, ec, maxiter) {
    var i, f, prec = StarJs.Kepler['DEFAULT_PRECISION'];

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    m = StarJs.Math.mod(m, StarJs.Math.PI2);

    /* Iterative solution of Kepler equation with Newthon's method.
     */
    //var e0 = 0; /* TODO: initial value depends on eccentricity */
    /* Gary R. Smith.  A simple, efficient starting value for the
     * iterative solution of Kepler's equation.  // Celestial
     * Mechanics and Dynamical Astronomy.  Volume 19, Number 2 (1979)
     * DOI: 10.1007/BF01796088.
     */
    var sinm = Math.sin(m);
    var e0 = m + ec * sinm / (1 - Math.sin(m + ec) + sinm);

    do {
        f = e0 - ec * Math.sin(e0) - m;
        e0 -= f / (1.0 - ec * Math.cos(e0));
    } while (maxiter-- > 0 && (f > prec));
    return (maxiter > 0) ? e0 : null;
};

/** Compute position of a body on elliptical orbit.
 * @param {number} gm Gravity constant.
 * @param {number} m Mean anomaly.
 * @param {number} a Semi-major axis.
 * @param {number} ec Eccentricity.
 */
StarJs.Kepler.elliptic = function (gm, m, a, ec) {
    var k = Math.sqrt(gm / a), e = StarJs.Kepler.eccAnomaly(m, ec);
    var cosE = Math.cos(e), sinE = Math.sin(e);
    var fac = Math.sqrt((1.0 - ec) * (1.0 + ec));
    return new StarJs.Vector.Vector3(a * (cosE - ec), a * fac * sinE, 0.0);
    // var rho = 1.0 - ec * cosE;
    // var vel = StarJs.Vector.Vector3(-k * sinE / rho, k * fac * cosE / rho, 0.0);
};

/** Compute position of a body on parabolic orbit.
 @param {number} gm Gravity constant.
 @param {number} t0   time of pericenter
 @param {number} t    time to calculate position for
 @param {number} q    pericenter distance
 @param {number} ec Eccentricity.
 @param {number=} maxiter Optional maximal number of iterations.
 */
StarJs.Kepler.parabolic = function (gm, t0, t, q, ec, maxiter) {
    function stumpff(e2, ret) {
        var eps = StarJs.Kepler['DEFAULT_PRECISION'];
        var n, add, c1, c2, c3;

        c1 = c2 = c3 = 0.0; 
        add = n = 1.0;

        do {
            c1 += add;  add /= (2.0*n);
            c2 += add;  add /= (2.0*n+1.0);
            c3 += add;  add *= -e2; 
            n += 1.0;
        }
        while (Math.abs(add) >= eps);

        ret.c1 = c1;
        ret.c2 = c2;
        ret.c3 = c3;
    }

    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }

    var e2 = 0, e20, fac, c = {}, k, tau, a, b, u, u2, r;
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    fac = 0.5 * ec;
    k   = Math.sqrt(gm/(q*(1+ec)));
    tau = Math.sqrt(gm)*(t-t0);
    do {
        e20 = e2;
        a = 1.5 * Math.sqrt(fac/(q*q*q))*tau;
        b = Math.pow(Math.sqrt(a*a + 1)+a, 1/3.0);
        u = b - 1.0/b;
        u2 = u*u;
        e2 = u2*(1-ec)/fac;
        stumpff(e2, c);
        fac = 3.0 * ec * c.c3;
    } while (Math.abs(e2-e20) >= prec);
    return new StarJs.Vector.Vector3(q*(1-u2*c.c2/fac),
                                     q*Math.sqrt((1+ec)/fac)*u*c.c1,
                                     0);
//     // res is position
//     r = q * (1+u2*c.c2*ec/fac);
//     var vel = new StarJs.Vector.Vector3(-k*res.y/r,
//                                         k*(res.x/r+ec),
//                                         0.0);
};

/** Compute hyperbolic anomaly.
 * @param {number} mh
 * @param {number} e
 * @param {number=} maxiter
 */
StarJs.Kepler.hypAnom = function (mh, e, maxiter) {
    if (maxiter === undefined) {
        maxiter = StarJs.Kepler['DEFAULT_ITERATIONS'];
    }
    var prec = StarJs.Kepler['DEFAULT_PRECISION'];

    var i = 0, h, f;

    h = Math.log(2.0*Math.abs(mh)/e+1.8);
    if (mh < 0.0) h = -h;

    do {
        f = e*StarJs.Math.sinh(h) - h - mh;
        h -= f/(e*StarJs.Math.cosh(h) - 1.0);
        ++i;
        if (i === maxiter) {
            // TODO: throw exception?
            return null;
        }
    } while (Math.abs(f) > prec*(1.0 + Math.abs(h+mh)));

    return h;
};

/** Compute position of a body on hyperbolic orbit.
 * @param {number} gm Gravity constant
 * @param {number} t0 Time of pericenter
 * @param {number} t Time
 * @param {number} a Semi-major axis
 * @param {number} e Eccentricity.
 */
StarJs.Kepler.hyperbolic = function (gm, t0, t, a, e) {
    var k, mh, h, ch, sh, rho, fac;

    a = Math.abs(a);
    k = Math.sqrt(gm/a);
    mh = k*(t-t0)/a;
    h = StarJs.Kepler.hypAnom(mh, e);
    ch = StarJs.Math.cosh(h);
    sh = StarJs.Math.sinh(h);
    fac = Math.sqrt((e+1)*(e-1)); // (e*e-1) ?
    rho = e*ch - 1;
    return new StarJs.Vector.Vector3(a*(e-ch), a*fac*sh, 0);
//     // ret is position
//     vel = new StarJs.Vector.Vector3(-k*sh/rho, k*fac*ch/rho, 0);
};

/** Calculate Keplerian orbital position.
    @param {number} gm   GM
    @param {number} t0   time of pericenter
    @param {number} t    time to calculate position for
    @param {number} q    pericenter distance
    @param {number} e    eccentricity
    @param {number} pqr  Gauss' vector matrix
 */
StarJs.Kepler.keplerPos = function (gm, t0, t, q, e, pqr) {
    var M0 = 0.1, EPS = 0.1, m, delta, tau, invax, r;

    delta = Math.abs(1-e);
    invax = delta / q;
    tau = Math.sqrt(gm)*(t-t0);
    m = tau * Math.sqrt(invax*invax*invax);

    if ((m < M0) && (delta < EPS)) {
        r = StarJs.Kepler.parabolic(gm, t0, t, q, e);
    } else if (e < 1.0) {
        r = StarJs.Kepler.elliptic(gm, m, 1.0/invax, e);
    } else {
        r = StarJs.Kepler.hyperbolic(gm, t0, t, 1.0/invax, e);
    }

    return pqr.apply(r);
};

/** Return Gauss vectors matrix for given elements.
    @param {number} omega longitude of the ascending node (radians)
    @param {number} i     inclination (radians)
    @param {number} w     argument of pericenter (radians)
 */
StarJs.Kepler.gaussVec = function (omega, i, w) {
    return StarJs.Vector.Matrix3.r_z(-omega).mult(StarJs.Vector.Matrix3.r_x(-i))
        .mult(StarJs.Vector.Matrix3.r_z(-w));
};
StarJs.Solar = {};

/**
 * Earth orbit's axial tilt.
 * @const
 * @type {number}
 */
StarJs.Solar.EPS = 23.43920111 * StarJs.Math.DEG2RAD;

/**
 * Return approximate position of Moon, quickly.
 * @param {number} mct Julian centuries from J2000. TODO
 * @return Object with fields ra and dec.
 */
StarJs.Solar.approxMoon = function (mct) {
    var l0, l, ls, f, d, dl, s, h, n, ml, mb, me;

    l0 = StarJs.Math.frac(0.606433 + 1336.855225 * mct);
    l  = StarJs.Math.PI2 * StarJs.Math.frac(0.374897 + 1325.552410 * mct);
    ls = StarJs.Math.PI2 * StarJs.Math.frac(0.993133 + 99.997361 * mct);
    d  = StarJs.Math.PI2 * StarJs.Math.frac(0.827361 + 1236.853086 * mct);
    f  = StarJs.Math.PI2 * StarJs.Math.frac(0.259086 + 1342.227825 * mct);

    dl = +22640 * Math.sin(l) - 4586 * Math.sin(l - 2 * d) + 2370 * Math.sin(2 * d) +
        769 * Math.sin(2 * l) - 668 * Math.sin(ls) - 412 * Math.sin(2 * f) +
        -212 * Math.sin(2 * l - 2 * d) - 206 * Math.sin(l + ls - 2 * d) +
        192 * Math.sin(l + 2 * d) - 165 * Math.sin(ls - 2 * d) +
        -125 * Math.sin(d) - 110 * Math.sin(l + ls) + 148 * Math.sin(l - ls) +
        -55 * Math.sin(2 * f - 2 * d);
    s = f + (dl + 412 * Math.sin(2 * f) + 541 * Math.sin(ls)) * StarJs.Math.ARCS2RAD;
    h = f - 2 * d;
    n = -526 * Math.sin(h) + 44 * Math.sin(l + h) - 31 * Math.sin(h - l) +
        -23 * Math.sin(ls + h) + 11 * Math.sin(h - ls) - 25 * Math.sin(f - 2 * l) +
        21 * Math.sin(f - l);
    ml = StarJs.Math.PI2 * StarJs.Math.frac(l0 + dl / 1296.0e3);
    mb = (18520.0 * Math.sin(s) + n) * StarJs.Math.ARCS2RAD;
    me = StarJs.Vector.Matrix3.r_x(-StarJs.Solar.EPS).apply((new StarJs.Vector.Polar3(ml, mb)).toVector3());
    me = new StarJs.Vector.Polar3(me);
    return {'ra': me['phi'], 'dec': me['theta']};

};

/**
 * Return approximate position of Moon, quickly.
 * @param {number} mct Julian centuries from J2000.
 * @return Object with fields ra and dec.
 */
StarJs.Solar.approxSun = function (mct) {
    var l, m, m2, me, se;
    m2 = StarJs.Math.frac(0.993133 + 99.997361 * mct);
    m = StarJs.Math.PI2 * m2;
    l = StarJs.Math.PI2 * StarJs.Math.frac(
        0.7859453  + m2 + (6893.0 * Math.sin(m) +
                           72.0 * Math.sin(2 * m) +
                           6191.2 * mct) / 1296.0e3);
    me = StarJs.Vector.Matrix3.r_x(-StarJs.Solar.EPS).apply((new StarJs.Vector.Polar3(l, 0)).toVector3());
    me = new StarJs.Vector.Polar3(me);
    return {'ra': me['phi'], 'dec': me['theta']};
};

/** @const */
StarJs.Solar.EVENTS = [{
    'body': 'moon',
    'name': 'day',
    'title': "Moon",
    'sinh0': Math.sin((+8.0 / 60.0) * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxMoon
}, {
    'body': 'sun',
    'name': 'day',
    'title': "Sun",
    'sinh0': Math.sin((-50.0 / 60.0) * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightC',
    'title': "Civil twilight",
    'sinh0': Math.sin(-6.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightN',
    'title': "Nautical twilight",
    'sinh0': Math.sin(-12.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}, {
    'body': 'sun',
    'name': 'twilightA',
    'title': "Astronomical twilight",
    'sinh0': Math.sin(-18.0 * StarJs.Math.DEG2RAD),
    'posFunc': StarJs.Solar.approxSun
}];

/** Calculate Sun and Moon events: rise, set and twilight.
 *
 * @param {number} start start time (MJD)
 * @param {number} end end day (MJD)
 * @param {number} lambda geographical longitude
 * @param {number} phi geographical latittude
 * @param {function(number)=} tz timezone offset function (converts UTC to local).
 *
 * Returns array of objects, each object describes particular day in form:
 *
 *  {
 *    moon: { 
 *      midnight: -0.4575,
 *      day: {
 *          rise: 10.2689,
 *          set: 21.0516
 *      }
 *    },
 *    sun:  { 
 *      midnight: -0.3877,
 *      day:       { rise: 6.3015, set: 20.56 },
 *      twilightA: { rise: 3.9416, set: 22.9039 },
 *      twilightN: { ... },
 *      twilightC: { ... }
 *    }
 *  }
 *
 * 'moon' and 'sun' objects have property 'midnight' that gives sinAlt
 * value of Moon and Sun at beginning of the day.  If rise and set
 * values absent, boolean 'alwaysAbove' field is set.
 */
StarJs.Solar.sunAndMoonEvents = function (start, end, lambda, phi, tz) {
    function sinAlt(approxBody, mjd, lambda, cphi, sphi) {
        var t, pos, tau;
        t = (mjd - 51544.5) / 36525.0;
        pos = approxBody(t);
        tau = StarJs.Time.gmst(mjd) + lambda - pos['ra'];
        return sphi * Math.sin(pos['dec']) +
            cphi * Math.cos(pos['dec']) * Math.cos(tau);
    }

    /** @const */
    var EVENTS = StarJs.Solar.EVENTS;

    var result = [], today, h, j, ptime, ctime, ntime, H = 1.0 / 24.0, posp = {}, pos0 = {}, posn = {}, cphi = Math.cos(phi), sphi = Math.sin(phi), name, evt, yp = {'sun': {}, 'moon': {}}, y0 = {'sun': {}, 'moon': {}}, yn = {'sun': {}, 'moon': {}}, interp = {};
    if (typeof tz === 'undefined') {
        tz = function (a) {
            return a;
        };
    }

    ptime = start;

    posp['moon'] = sinAlt(StarJs.Solar.approxMoon, ptime, lambda, cphi, sphi);
    posp['sun']  = sinAlt(StarJs.Solar.approxSun,  ptime, lambda, cphi, sphi);

    for (j = 0; j < EVENTS.length; j += 1) {
        evt = EVENTS[j];
        name = evt['name'];

        yp[evt['body']][name] = posp[evt['body']] - evt['sinh0'];
    }

    while (ptime < end) {
        today = {
            'moon': {
                'midnight': posp['moon']
            },
            'sun': {
                'midnight': posp['sun']
            }
        };
        for (h = 1; h < 24; h += 2) {
            ctime = ptime + H;
            ntime = ctime + H; // ntime = ctime + 2 * H;

            // Calc new positions...
            pos0['moon'] = sinAlt(StarJs.Solar.approxMoon, ctime, lambda,
                               cphi, sphi);
            pos0['sun'] = sinAlt(StarJs.Solar.approxSun, ctime, lambda,
                              cphi, sphi);
            posn['moon'] = sinAlt(StarJs.Solar.approxMoon, ntime, lambda,
                               cphi, sphi);
            posn['sun'] = sinAlt(StarJs.Solar.approxSun, ntime, lambda,
                              cphi, sphi);

            for (j = 0; j < EVENTS.length; j += 1) {
                evt = EVENTS[j];
                name = evt['name'];
                today[evt['body']][name] = today[evt['body']][name] || {};

                y0[evt['body']][name] = pos0[evt['body']] - evt['sinh0'];
                yn[evt['body']][name] = posn[evt['body']] - evt['sinh0'];

                interp = StarJs.Math.quadInterpolation(yp[evt['body']][name], y0[evt['body']][name], yn[evt['body']][name]);

                switch (interp['roots'].length) {
                case 0:
                    // No roots
                    break;
                case 1:
                    // Single root
                    if (yp[evt['body']][name] < 0.0) {
                        today[evt['body']][name]['rise'] = h + interp['roots'][0];
                    } else {
                        today[evt['body']][name]['set']  = h + interp['roots'][0];
                    }
                    break;
                case 2:
                    // Two roots
                    if (interp['ye'] < 0.0) {
                        today[evt['body']][name]['rise'] = h + interp['roots'][1];
                        today[evt['body']][name]['set']  = h + interp['roots'][0];
                    } else {
                        today[evt['body']][name]['rise'] = h + interp['roots'][0];
                        today[evt['body']][name]['set']  = h + interp['roots'][1];
                    }
                    break;
                }
            }

            // Next interval
            yp = yn;
            yn = {'moon': {}, 'sun': {}};

            ptime = ntime;
        }

        for (j = 0; j < EVENTS.length; j += 1) {
            evt = EVENTS[j];
            name = evt['name'];

            if (!today[evt['body']][name]['set'] && !today[evt['body']][name]['rise']) {
                today[evt['body']][name]['alwaysAbove'] = (pos0[evt['body']] > evt['sinh0']);
            }
        }
        // Next day
        ptime = (start += 1.0);

        result.push(today);
    }

    return result;
};

/**********************************************************************
 *
 * Solar sytem plantes.
 *
 */

/** @constructor */
StarJs.Solar.Sun = function () {
    this['name'] = 'Sun';
};

/** @const */
StarJs.Solar.Sun.$POS = new StarJs.Vector.Vector3(0, 0, 0);

/** @const */
StarJs.Solar.GM = StarJs.Math.sqr(0.01720209895);

/**
 * Cooordinates of Sun are constant.
 */
StarJs.Solar.Sun.prototype['keplerCoord'] = function (t) {
    return StarJs.Solar.Sun.$POS;
};

/** @constructor */
StarJs.Solar.Body = function (name, params) {
    this['name'] = name;
    this['a']  = params['a'];
    this['ec'] = params['e'];
    this['m0'] = params['M0'];
    this['n']  = params['n'];
    this['omega'] = params['O'];
    this['i']  = params['i'];
    this['w']  = params['w'];
    this['t0'] = params['T0'];
    // ...
};

StarJs.Solar.Body.prototype['keplerCoord'] = function (t) {
    var P = 1.3970;
    var m = StarJs.Math.DEG2RAD * (this['m0'] + this['n'] * (t - this['t0']));
    var r = StarJs.Kepler.elliptic(StarJs.Solar.GM, m, this['a'], this['ec']);
    var pqr = StarJs.Kepler.gaussVec(
        StarJs.Math.DEG2RAD * (this['omega'] + P * t),
        StarJs.Math.DEG2RAD * this['i'],
        StarJs.Math.DEG2RAD * (this['w'] - this['omega']));
    return pqr.apply(r);
};

/** @const */
StarJs.Solar.BODIES = {
    'Sun':     new StarJs.Solar.Sun(),
    'Mercury': new StarJs.Solar.Body('Mercury', {
      'a' :  0.387099, 'e' : 0.205634, 'M0' : 174.7947, 'n'  : 149472.6738,
      'O' :  48.331,   'i' : 7.0048,   'w'  :  77.4552, 'T0' : 0.0
    }),
    'Venus':   new StarJs.Solar.Body('Venus', {
      'a' :  0.723332, 'e' : 0.006773, 'M0' :  50.4071, 'n'  : 58517.8149,
      'O' :  76.680,   'i' : 3.3946,   'w'  : 131.5718, 'T0' : 0.0
    }),
    'Earth':   new StarJs.Solar.Body('Earth', {
      'a' :  1.000000, 'e' : 0.016709, 'M0' : 357.5256, 'n'  : 35999.3720,
      'O' : 174.876,   'i' : 0.0000,   'w'  : 102.9400, 'T0' : 0.0
    }),
    'Mars':   new StarJs.Solar.Body('Mars', {
      'a' :  1.523692, 'e' : 0.093405, 'M0' :  19.3879, 'n'  : 19140.3023,
      'O' :  49.557,   'i' : 1.8496,   'w'  : 336.0590, 'T0' : 0.0
    }),
    'Jupiter':   new StarJs.Solar.Body('Jupiter', {
      'a' :  5.204267, 'e' : 0.048775, 'M0' :  18.8185, 'n'  : 3033.6272,
      'O' : 100.4908,  'i' : 1.3046,   'w'  :  15.5576, 'T0' : 0.0
    }),
    'Saturn':   new StarJs.Solar.Body('Saturn', {
      'a' :  9.582018, 'e' : 0.055723, 'M0' : 320.3477, 'n'  : 1213.8664,
      'O' : 113.6427,  'i' : 2.4852,   'w'  :  89.6567, 'T0' : 0.0
    }),
    'Uranus':   new StarJs.Solar.Body('Uranus', {
      'a' : 19.229412, 'e' : 0.044406, 'M0' : 142.9559, 'n'  : 426.9282,
      'O' :  73.9893,  'i' : 0.7726,   'w'  : 170.5310, 'T0' : 0.0
    }),
    'Neptune':   new StarJs.Solar.Body('Neptune', {
      'a' : 30.103658, 'e' : 0.011214, 'M0' : 267.7649, 'n'  : 217.9599,
      'O' : 131.7942,  'i' : 1.7680,   'w'  :  37.4435, 'T0' : 0.0
    })
};
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
    (new StarMap.ConstellationBoundaries(window['CON_BOUND_18'])).draw(ctx, this['proj']);

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
