/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <numeric>
#include "helper_functions.h"

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // Initialize noise generators
    gen = default_random_engine();
    n_x = normal_distribution<double>(0, std[0]);
    n_y = normal_distribution<double>(0, std[1]);
    n_theta = normal_distribution<double>(0, std[2]);

    // Initialize all particles
    num_particles = 250000;
    particles.resize(num_particles);
    weights.resize(num_particles);
    for (int i = 0; i < num_particles; i++) {
        particles[i].id = i;
        particles[i].x = x + n_x(gen);
        particles[i].y = y + n_y(gen);
        particles[i].theta = theta + n_theta(gen);
        particles[i].weight = 1;
        weights[i] = 1;
    }

    is_initialized = true;
    cout << "Finished Init!!!" << endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    // Reset noise generators to passed in noise tolerances.
    updateNoiseGenerators(std_pos);

    // Apply motion model to each particle.
    double theta_orig;
    for (int p = 0; p < particles.size(); p++) {
        if (yaw_rate < 0.001) {
            particles[p].x += velocity * delta_t * cos(particles[p].theta) + n_x(gen);
            particles[p].y += velocity * delta_t * sin(particles[p].theta) + n_y(gen);
        } else {
            theta_orig = particles[p].theta;
            particles[p].theta += yaw_rate * delta_t + n_theta(gen);
            particles[p].x += (velocity / yaw_rate) * (sin(particles[p].theta) - sin(theta_orig)) + n_x(gen);
            particles[p].y += (velocity / yaw_rate) * (cos(theta_orig) - cos(particles[p].theta)) + n_y(gen);
        };
    }


}

vector<int> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    // Vector to hold distance between observation and nearest landmark.
    vector<int> nearest_landmark_idx;
    nearest_landmark_idx.resize(observations.size());
    double distance, min_distance;

    for (int i=0; i < observations.size(); i++) {
        nearest_landmark_idx[i] = 0;
        min_distance = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
        for (int j=1; j < predicted.size(); j++) {
            distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if (distance < min_distance) {
                nearest_landmark_idx[i] = j;
                min_distance = distance;
            }
        }
    }
    return nearest_landmark_idx;

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
    //   for the fact that the map's y-axis actually points downwards.)
    //   http://planning.cs.uiuc.edu/node99.html

    // Vector to hold transformed observations
    vector<LandmarkObs> t_obs;
    t_obs.resize(observations.size());

    // Loop over particles
    for (int i = 0; i < particles.size(); i++) {

        // Transform observations to map-coordinates assuming current particle pose
        for (int j = 0; j < t_obs.size(); j++) {
            t_obs[j].x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) +
                         particles[i].x;
            t_obs[j].y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) +
                         particles[i].y;
        }

        // Find all landmarks within sensor_range of the current particle
        vector<LandmarkObs> nearby_landmarks;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double landmark_dist = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f,
                                        map_landmarks.landmark_list[j].y_f);
            //cout << landmark_dist << endl;
            if (landmark_dist <= sensor_range) {
                LandmarkObs nearby_landmark;
                nearby_landmark.x = map_landmarks.landmark_list[j].x_f;
                nearby_landmark.y = map_landmarks.landmark_list[j].y_f;
                nearby_landmarks.push_back(nearby_landmark);
                //cout << "Pushed!!" << endl;
            }

        }


        // For each transformed observation, find distance to nearest landmark within sensor range.
        vector<int> nearest_landmark_idx = dataAssociation(nearby_landmarks, t_obs);


        // Compute particle weight by evaluating multivariate gaussian using observations and their nearest landmarks.
        particles[i].weight = 1;
        for (int j=0; j < t_obs.size(); j++) {
            particles[i].weight *= evaluteMVN(t_obs[j].x, t_obs[j].y,
                                              nearby_landmarks[nearest_landmark_idx[j]].x,
                                              nearby_landmarks[nearest_landmark_idx[j]].y,
                                              std_landmark[0], std_landmark[1]);
        }

        // Update vector of all particle weights
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::discrete_distribution<int> sampler (weights.begin(), weights.end());
    vector<Particle> sampled_particles;
    sampled_particles.resize(particles.size());
    for (int i = 0; i < num_particles; i++) {
        sampled_particles[i] = particles[sampler(gen)];
    }
    particles.assign(sampled_particles.begin(), sampled_particles.end());
}

void ParticleFilter::write(std::string filename) {
    // You don't need to modify this file.
    std::ofstream dataFile;
    dataFile.open(filename, std::ios::app);
    for (int i = 0; i < num_particles; ++i) {
        dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    dataFile.close();
}

void ParticleFilter::updateNoiseGenerators(double *std) {
    n_x.param(normal_distribution<double>::param_type(0, std[0]));
    n_y.param(normal_distribution<double>::param_type(0, std[1]));
    n_theta.param(normal_distribution<double>::param_type(0, std[2]));
}

double ParticleFilter::evaluteMVN(double x, double y, double mu_x, double mu_y, double sigma_x, double sigma_y) {
    const double pi = 3.1415926535897;
    double weight = (1/(2*pi*sigma_x*sigma_y)) * exp(-0.5*( pow(x-mu_x,2)/pow(sigma_x,2) + pow(y-mu_y,2)/pow(sigma_y,2)));
    return weight;
}
