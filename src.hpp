#ifndef PPCA_SRC_HPP
#define PPCA_SRC_HPP
#include "math.h"
class Monitor;

class Controller {

public:
    Controller(const Vec &_pos_tar, double _v_max, double _r, int _id, Monitor *_monitor) {
        pos_tar = _pos_tar;
        v_max = _v_max;
        r = _r;
        id = _id;
        monitor = _monitor;
    }

    void set_pos_cur(const Vec &_pos_cur) {
        pos_cur = _pos_cur;
    }

    void set_v_cur(const Vec &_v_cur) {
        v_cur = _v_cur;
    }

private:
    int id;
    Vec pos_tar;
    Vec pos_cur;
    Vec v_cur;
    double v_max, r;
    Monitor *monitor;

    // Helper: check if a planned velocity would collide with any other robot during next interval
    bool will_collide_with_any(const Vec &v_plan) const {
        int n = monitor->get_robot_number();
        for (int j = 0; j < n; ++j) {
            if (j == id) continue;
            // Yield to lower-id robots only
            if (j < id && will_collide_with(j, v_plan)) return true;
        }
        return false;
    }

    bool predict_collision(const Vec &v_plan, int &other_id_out) const {
        int n = monitor->get_robot_number();
        for (int j = 0; j < n; ++j) {
            if (j == id) continue;
            if (j < id && will_collide_with(j, v_plan)) { other_id_out = j; return true; }
        }
        return false;
    }

    bool will_collide_with(int other_id, const Vec &v_plan) const {
        Vec other_pos = monitor->get_pos_cur(other_id);
        Vec other_v = monitor->get_v_cur(other_id);
        double other_r = monitor->get_r(other_id);

        Vec delta_pos = pos_cur - other_pos;
        Vec delta_v = v_plan - other_v;
        double delta_v_norm = delta_v.norm();
        double project = delta_pos.dot(delta_v);
        if (project >= 0) return false;
        project /= -delta_v_norm;
        double delta_r = r + other_r;
        double min_dis_sqr;
        if (project < delta_v_norm * TIME_INTERVAL) {
            min_dis_sqr = delta_pos.norm_sqr() - project * project;
        } else {
            min_dis_sqr = (delta_pos + delta_v * TIME_INTERVAL).norm_sqr();
        }
        return min_dis_sqr <= delta_r * delta_r - EPSILON;
    }

    Vec clamp_speed(const Vec &v) const {
        double speed = v.norm();
        double vmax = std::max(0.0, v_max - 1e-6);
        if (speed <= vmax) return v;
        return v.normalize() * vmax;
    }

public:

    Vec get_v_next() {
        // Base desired direction towards target
        Vec to_tar = pos_tar - pos_cur;
        double dist = to_tar.norm();
        if (dist <= EPSILON) return Vec();

        // Prefer a speed that reaches target in one step but capped by v_max
        double desired_speed = std::min(v_max, dist / TIME_INTERVAL);

        // Keep full speed; do not globally slow down on warnings to avoid stalling

        Vec dir = to_tar.normalize();
        Vec v_des = dir * desired_speed;
        // If previous step had any warning globally, bias a bit sideways to avoid repeating conflicts
        if (monitor->get_warning()) {
            Vec tilt(-dir.y, dir.x);
            double side = (id % 2 == 0) ? 1.0 : -1.0;
            v_des = clamp_speed(v_des + tilt * (desired_speed * 0.2 * side));
        }

        // Add repulsion from nearby robots to break symmetry and maintain separation
        Vec rep(0, 0);
        int n = monitor->get_robot_number();
        for (int j = 0; j < n; ++j) {
            if (j == id) continue;
            Vec o_pos = monitor->get_pos_cur(j);
            double o_r = monitor->get_r(j);
            Vec delta = pos_cur - o_pos;
            double d = delta.norm();
            double sep = r + o_r + std::max(0.5, v_max * TIME_INTERVAL * 1.5);
            if (d < sep) {
                double t = (sep - d) / sep; // in (0,1]
                double bias = (j < id) ? 1.4 : 0.8; // yield more to lower-id robots
                rep += delta.normalize() * (desired_speed * bias * t);
            }
        }

        Vec v_plan = clamp_speed(v_des + rep);

        // Fast path: if no predicted collision, use it
        if (!will_collide_with_any(v_plan)) return v_plan;

        // Try to reduce speed gradually
        double factor = 0.8;
        for (int tries = 0; tries < 6; ++tries) {
            Vec v_try = clamp_speed(dir * (desired_speed * factor));
            if (!will_collide_with_any(v_try)) return v_try;
            factor *= 0.7;
        }

        // Try slight sidestep by rotating direction with small angle depending on id parity
        double sign = (id % 2 == 0) ? 1.0 : -1.0;
        double angles[10] = {0.15, 0.3, 0.45, 0.6, 0.9, -0.15, -0.3, -0.45, -0.6, -0.9};
        for (double ang : angles) {
            Vec d = dir.rotate(sign * ang);
            Vec v_try = clamp_speed(d * (desired_speed * 0.6));
            if (!will_collide_with_any(v_try)) return v_try;
        }

        // If still colliding with a specific nearest robot, try moving tangentially relative to that robot
        int other_id = -1;
        if (predict_collision(v_plan, other_id) && other_id >= 0) {
            Vec o_pos = monitor->get_pos_cur(other_id);
            Vec radial = (pos_cur - o_pos).normalize();
            Vec tangent = radial.rotate(PI / 2);
            Vec v_try = clamp_speed((dir * (desired_speed * 0.5)) + (tangent * (desired_speed * 0.6)));
            if (!will_collide_with_any(v_try)) return v_try;
            v_try = clamp_speed((dir * (desired_speed * 0.5)) - (tangent * (desired_speed * 0.6)));
            if (!will_collide_with_any(v_try)) return v_try;
            // fallback to pure tangent small move
            v_try = clamp_speed(tangent * (desired_speed * 0.3));
            if (!will_collide_with_any(v_try)) return v_try;
            v_try = clamp_speed(tangent * (-desired_speed * 0.3));
            if (!will_collide_with_any(v_try)) return v_try;
        }

        // As a near-last resort, drift slightly away from neighbors to break deadlock
        Vec repel(0, 0);
        int nn = monitor->get_robot_number();
        for (int j = 0; j < nn; ++j) {
            if (j == id) continue;
            Vec o_pos = monitor->get_pos_cur(j);
            repel += (pos_cur - o_pos);
        }
        if (repel.norm() > 1e-6) {
            Vec v_try = clamp_speed(repel.normalize() * std::max(0.1 * v_max, desired_speed * 0.3));
            if (!will_collide_with_any(v_try)) return v_try;
        }

        // As a last resort, stop to avoid collision
        return Vec();
    }
};


#endif //PPCA_SRC_HPP
