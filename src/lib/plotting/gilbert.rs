//! These functions are a port of those in https://github.com/jakubcerveny/gilbert, used under
//! the following license:
//!
//! BSD 2-Clause License
//!
//! Copyright (c) 2018, Jakub Červený
//! All rights reserved.
//!
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions are met:
//!
//! * Redistributions of source code must retain the above copyright notice, this
//!   list of conditions and the following disclaimer.
//!
//! * Redistributions in binary form must reproduce the above copyright notice,
//!   this list of conditions and the following disclaimer in the documentation
//!   and/or other materials provided with the distribution.
//!
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/// Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
/// 2D rectangular grids. Takes a position along the gilbert curve and returns
/// its 2D (x,y) coordinate.
pub fn gilbert_d2xy(idx: u32, width: u32, height: u32) -> (u32, u32) {
    let (x, y) = if width >= height {
        gilbert_d2xy_r(idx as i32, 0, 0, 0, width as i32, 0, 0, height as i32)
    } else {
        gilbert_d2xy_r(idx as i32, 0, 0, 0, 0, height as i32, width as i32, 0)
    };
    (x as u32, y as u32)
}

fn sgn(x: i32) -> i32 {
    match x.cmp(&0) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

fn gilbert_d2xy_r(
    dst_idx: i32,
    mut cur_idx: i32,
    x: i32,
    y: i32,
    ax: i32,
    ay: i32,
    bx: i32,
    by: i32,
) -> (i32, i32) {
    let w = (ax + ay).abs();
    let h = (bx + by).abs();

    let (dax, day) = (sgn(ax), sgn(ay)); // unit major direction
    let (dbx, dby) = (sgn(bx), sgn(by)); // unit orthogonal direction

    let _dx = dax + dbx;
    let _dy = day + dby;
    let di = dst_idx - cur_idx;

    if h == 1 {
        return (x + dax * di, y + day * di);
    }
    if w == 1 {
        return (x + dbx * di, y + dby * di);
    }

    let (mut ax2, mut ay2) = (ax / 2, ay / 2);
    let (mut bx2, mut by2) = (bx / 2, by / 2);

    let w2 = (ax2 + ay2).abs();
    let h2 = (bx2 + by2).abs();

    if 2 * w > 3 * h {
        if w2 % 2 != 0 && w > 2 {
            // prefer even steps
            (ax2, ay2) = (ax2 + dax, ay2 + day)
        }

        // long case: split in two parts only
        let nxt_idx = cur_idx + ((ax2 + ay2) * (bx + by)).abs();
        if cur_idx <= dst_idx && dst_idx < nxt_idx {
            return gilbert_d2xy_r(dst_idx, cur_idx, x, y, ax2, ay2, bx, by);
        }
        cur_idx = nxt_idx;

        return gilbert_d2xy_r(
            dst_idx,
            cur_idx,
            x + ax2,
            y + ay2,
            ax - ax2,
            ay - ay2,
            bx,
            by,
        );
    }

    if h2 % 2 != 0 && h > 2 {
        // prefer even steps
        (bx2, by2) = (bx2 + dbx, by2 + dby);
    }

    // standard case: one step up, one long horizontal, one step down
    let nxt_idx = cur_idx + ((bx2 + by2) * (ax2 + ay2)).abs();
    if cur_idx <= dst_idx && dst_idx < nxt_idx {
        return gilbert_d2xy_r(dst_idx, cur_idx, x, y, bx2, by2, ax2, ay2);
    }
    cur_idx = nxt_idx;

    let nxt_idx = cur_idx + ((ax + ay) * ((bx - bx2) + (by - by2))).abs();
    if cur_idx <= dst_idx && dst_idx < nxt_idx {
        return gilbert_d2xy_r(
            dst_idx,
            cur_idx,
            x + bx2,
            y + by2,
            ax,
            ay,
            bx - bx2,
            by - by2,
        );
    }
    cur_idx = nxt_idx;

    gilbert_d2xy_r(
        dst_idx,
        cur_idx,
        x + (ax - dax) + (bx2 - dbx),
        y + (ay - day) + (by2 - dby),
        -bx2,
        -by2,
        -(ax - ax2),
        -(ay - ay2),
    )
}
