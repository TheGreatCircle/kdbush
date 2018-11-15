(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
typeof define === 'function' && define.amd ? define(factory) :
(global.KDBush = factory());
}(this, (function () { 'use strict';

function sortKD(ids, coords, nodeSize, left, right, axis) {
    if (right - left <= nodeSize) { return; }

    var m = (left + right) >> 1; // middle index

    // sort ids and coords around the middle index so that the halves lie
    // either left/right or top/bottom correspondingly (taking turns)
    select(ids, coords, m, left, right, axis);

    // recursively kd-sort first half and second half on the opposite axis
    sortKD(ids, coords, nodeSize, left, m - 1, 1 - axis);
    sortKD(ids, coords, nodeSize, m + 1, right, 1 - axis);
}

// custom Floyd-Rivest selection algorithm: sort ids and coords so that
// [left..k-1] items are smaller than k-th item (on either x or y axis)
function select(ids, coords, k, left, right, axis) {

    while (right > left) {
        if (right - left > 600) {
            var n = right - left + 1;
            var m = k - left + 1;
            var z = Math.log(n);
            var s = 0.5 * Math.exp(2 * z / 3);
            var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            select(ids, coords, k, newLeft, newRight, axis);
        }

        var t = coords[2 * k + axis];
        var i = left;
        var j = right;

        swapItem(ids, coords, left, k);
        if (coords[2 * right + axis] > t) { swapItem(ids, coords, left, right); }

        while (i < j) {
            swapItem(ids, coords, i, j);
            i++;
            j--;
            while (coords[2 * i + axis] < t) { i++; }
            while (coords[2 * j + axis] > t) { j--; }
        }

        if (coords[2 * left + axis] === t) { swapItem(ids, coords, left, j); }
        else {
            j++;
            swapItem(ids, coords, j, right);
        }

        if (j <= k) { left = j + 1; }
        if (k <= j) { right = j - 1; }
    }
}

function swapItem(ids, coords, i, j) {
    swap(ids, i, j);
    swap(coords, 2 * i, 2 * j);
    swap(coords, 2 * i + 1, 2 * j + 1);
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function range(ids, coords, minX, minY, maxX, maxY, nodeSize) {
    var stack = [0, ids.length - 1, 0];
    var result = [];

    // recursively search for items in range in the kd-sorted arrays
    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        // if we reached "tree node", search linearly
        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                var x = coords[2 * i];
                var y = coords[2 * i + 1];
                if (x >= minX && x <= maxX && y >= minY && y <= maxY) { result.push(ids[i]); }
            }
            continue;
        }

        // otherwise find the middle index
        var m = (left + right) >> 1;

        // include the middle item if it's in range
        var x$1 = coords[2 * m];
        var y$1 = coords[2 * m + 1];
        if (x$1 >= minX && x$1 <= maxX && y$1 >= minY && y$1 <= maxY) { result.push(ids[m]); }

        // queue search in halves that intersect the query
        if (axis === 0 ? minX <= x$1 : minY <= y$1) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(1 - axis);
        }
        if (axis === 0 ? maxX >= x$1 : maxY >= y$1) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(1 - axis);
        }
    }

    return result;
}

function within(ids, coords, qx, qy, r, nodeSize) {
    var stack = [0, ids.length - 1, 0];
    var result = [];
    var r2 = r * r;

    // recursively search for items within radius in the kd-sorted arrays
    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        // if we reached "tree node", search linearly
        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                if (sqDist(coords[2 * i], coords[2 * i + 1], qx, qy) <= r2) { result.push(ids[i]); }
            }
            continue;
        }

        // otherwise find the middle index
        var m = (left + right) >> 1;

        // include the middle item if it's in range
        var x = coords[2 * m];
        var y = coords[2 * m + 1];
        if (sqDist(x, y, qx, qy) <= r2) { result.push(ids[m]); }

        // queue search in halves that intersect the query
        if (axis === 0 ? qx - r <= x : qy - r <= y) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(1 - axis);
        }
        if (axis === 0 ? qx + r >= x : qy + r >= y) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(1 - axis);
        }
    }

    return result;
}

function sqDist(ax, ay, bx, by) {
    var dx = ax - bx;
    var dy = ay - by;
    return dx * dx + dy * dy;
}

var defaultGetX = function (p) { return p[0]; };
var defaultGetY = function (p) { return p[1]; };

var KDBush = function KDBush(
  points,
  ref
) {
  var getX = ref.getX; if ( getX === void 0 ) getX = defaultGetX;
  var getY = ref.getY; if ( getY === void 0 ) getY = defaultGetY;
  var nodeSize = ref.nodeSize; if ( nodeSize === void 0 ) nodeSize = 64;
  var ArrayType = ref.ArrayType; if ( ArrayType === void 0 ) ArrayType = Float64Array;
  var ids = ref.ids;
  var coords = ref.coords;

  this.nodeSize = 64;
  this.points = Float64Array;

  var IndexArrayType = points.length < 65536 ? Uint16Array : Uint32Array;

  if (ids && coords) {
    this.ids = ids;
    this.coords = coords;
  } else {
    // store indices to the input array and coordinates in separate typed arrays
    var ids$1 = (this.ids = new IndexArrayType(points.length));
    var coords$1 = (this.coords = new ArrayType(points.length * 2));

    for (var i = 0; i < points.length; i++) {
      ids$1[i] = i;
      coords$1[2 * i] = getX(points[i]);
      coords$1[2 * i + 1] = getY(points[i]);
    }

    // kd-sort both arrays for efficient search (see comments in sort.js)
    sortKD(ids$1, coords$1, nodeSize, 0, ids$1.length - 1, 0);
  }
};

KDBush.prototype.range = function range$1 (minX, minY, maxX, maxY) {
  return range(this.ids, this.coords, minX, minY, maxX, maxY, this.nodeSize)
};

KDBush.prototype.within = function within$1 (x, y, r) {
  return within(this.ids, this.coords, x, y, r, this.nodeSize)
};

return KDBush;

})));
