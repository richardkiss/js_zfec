var assert = require("assert");

var fec = require("./fec.js").fec;

function string_to_uint8(s) {
    var a = new Uint8Array(s.length);
    var i;
    for (i=0;i<s.length;i++) {
        a[i] = s.charCodeAt(i);
    }
    return a;
}

function uint8_to_string(ui) {
    var s, i;
    s = '';
    for (i=0;i<ui.length;i++) {
        s += String.fromCharCode(ui[i]);
    }
    return s;
}

function test_encode_decode(s) {
    var f = fec(3, 10);
    var encoded = f.encode(string_to_uint8(s), [0,1,2,3,4,5]);
    var i;
    var decoded;
    var tuples = [[0,1,2], [0,1,4], [0,1,3], [0,2,4], [0,1,5], [3,4,5], [2,3,4], [2,3,5], [2,4,5]];
    for (i=0;i<tuples.length;i++) {
        var tuple = tuples[i];
        decoded = f.decode([encoded[tuple[0]], encoded[tuple[1]], encoded[tuple[2]]], tuple);
        var s1 = uint8_to_string(decoded);
        assert.equal(s, s1)
    }
}

test_encode_decode("123456789");
test_encode_decode("0123456789");
test_encode_decode("0123456789abcdefghijklmnop");
test_encode_decode("Please kill me");
test_encode_decode("The quick brown fox jumps over the lazy dogs.");

var s = "abcdefghijklmnopqrstuvwxyz.!^*4&Y#&^(*&$(*&$(*%&@*($&%*($&%)))))";

while (s.length<163840) {
    test_encode_decode(s);
    s = s + s;
}
