/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <iostream>
#include "helper_functions.h"
#include "particle_filter.h"

using namespace std;

// Random generator or sampling.
default_random_engine gen;


/**
 * Evaluates 2D Multivariate Gaussian Distribution (MVN).
 * @param x Pointer to array corresponding to 2D point at which to evaluate MVN.
 * @param mu Pointer to array corresponding to 2D mean of MVN.
 * @param sigma Pointer to array corresponding to diagonal covariance matrix of MVN.
 * @return double corresponding to evaluationg MVN at argument x.
 */
double evaluateMVN(double *x, double *mu, double *sigma) {
    double mvn_value = 1.0;
    for (int i = 0; i < 2; i++) {
        mvn_value *= (1.0 / sigma[i]) * exp(-0.5 * pow(x[i] - mu[i], 2) / pow(sigma[i], 2));
    }
    mvn_value *= 1.0 / (2.0 * M_PI);
    return mvn_value;
}

/**
 * nearestLandmark Finds index of landmark nearest each observation.
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 * @return Vector index of landmark in predicted closest to each observation in observations.
 */
vector<int> nearestLandmark(std::vector<LandmarkObs> predicted,
                            std::vector<LandmarkObs> &observations) {

    // Vector to hold index of landmark in predicted that is closest to each observation.
    vector<int> nearest_landmark_idx;
    nearest_landmark_idx.resize(observations.size());

    // Variables to use while searching for landmark closest to each observation.
    double distance, min_distance;

    // For each observation, find landmark that is closest and save its index
    // in nearest_landmark_idx;
    for (int i = 0; i < observations.size(); i++) {
        nearest_landmark_idx[i] = 0;
        min_distance = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
        for (int j = 1; j < predicted.size(); j++) {
            distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if (distance < min_distance) {
                nearest_landmark_idx[i] = j;
                min_distance = distance;
            }
        }
    }
    return nearest_landmark_idx;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // Initialize noise generators for x, y, and theta.
    std::normal_distribution<double> n_x (0, std[0]);
    std::normal_distribution<double> n_y (0, std[1]);
    std::normal_distribution<double> n_theta (0, std[2]);

    // Initialize 1000 particles to GPS position with uncertainty and set all particle weights = 1
    num_particles = 1000;
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

    // Declare filter as being initialized.
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // Reset noise generators to supplied in noise sigma in argument std_pos.
    std::normal_distribution<double> n_x (0, std_pos[0]);
    std::normal_distribution<double> n_y (0, std_pos[1]);
    std::normal_distribution<double> n_theta (0, std_pos[2]);

    // Apply motion model to each particle.
    double theta_orig;
    for (int p = 0; p < particles.size(); p++) {
        // Update equations for yaw_rate > 0.
        if (fabs(yaw_rate) > 0.001) {
            theta_orig = particles[p].theta;
            particles[p].theta += yaw_rate * delta_t + n_theta(gen);
            particles[p].x += (velocity / yaw_rate) * (sin(particles[p].theta) - sin(theta_orig)) + n_x(gen);
            particles[p].y += (velocity / yaw_rate) * (cos(theta_orig) - cos(particles[p].theta)) + n_y(gen);
        }
        // Update equations for yaw_rate approximately 0.
        else {
            particles[p].x += velocity * delta_t * cos(particles[p].theta) + n_x(gen);
            particles[p].y += velocity * delta_t * sin(particles[p].theta) + n_y(gen);
            particles[p].theta += n_theta(gen);
        }
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.

    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    // NOTE: I used the method nearestLandmark (declared above) to perform data association.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // Vector to hold transformed observations
    vector<LandmarkObs> t_obs;
    t_obs.resize(observations.size());

    // Arrays to define 2D point and 2D mean for Multivariate Guassian (MVN) that will be used
    // to weight particles.
    double x[2];
    double mu[2];

    // Loop over each particle.
    for (int i = 0; i < particles.size(); i++) {

        // Transform observations to map-coordinates assuming current particle pose
        for (int j = 0; j < t_obs.size(); j++) {
            t_obs[j].x = observations[j].x * cos(particles[i].theta)  - observations[j].y * sin(particles[i].theta) +
                         particles[i].x;
            t_obs[j].y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) +
                         particles[i].y;
        }

        // Find all landmarks within sensor_range of the current particle
        vector<LandmarkObs> nearby_landmarks;
        double landmark_dist;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            landmark_dist = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f,
                                 map_landmarks.landmark_list[j].y_f);
            if (landmark_dist < sensor_range) {
                LandmarkObs nearby_landmark;
                nearby_landmark.x = map_landmarks.landmark_list[j].x_f;
                nearby_landmark.y = map_landmarks.landmark_list[j].y_f;
                nearby_landmarks.push_back(nearby_landmark);
            }

        }

        // If no landmarks are within sensor_range of current particle, then set current particle weight=0
        if (nearby_landmarks.size() > 0) {

            // For each transformed observation, find distance to nearest landmark within sensor range.
            vector<int> nearest_landmark_idx = nearestLandmark(nearby_landmarks, t_obs);

            // Compute particle weight by evaluating multivariate gaussian using each transformed observation
            // and the closest landmark.
            particles[i].weight = 1;
            for (int j = 0; j < t_obs.size(); j++) {
                x[0] = t_obs[j].x;
                x[1] = t_obs[j].y;
                mu[0] = nearby_landmarks[nearest_landmark_idx[j]].x;
                mu[1] = nearby_landmarks[nearest_landmark_idx[j]].y;
                particles[i].weight *= evaluateMVN(x, mu, std_landmark);;
            }
            weights[i] = particles[i].weight;
        } else {
            particles[i].weight = 0.0;
            weights[i] = particles[i].weight;
        }
    }
}

void ParticleFilter::resample() {
    // Discrete distribution from current particle weights.
    std::discrete_distribution<int> sampler(weights.begin(), weights.end());
    // Vector to hold particles that are sampled.
    vector<Particle> sampled_particles;
    sampled_particles.resize(particles.size());
    for (int i = 0; i < num_particles; i++) {
        sampled_particles[i] = particles[sampler(gen)];
    }
    // Set particles to newly sampled particles.
    particles.assign(sampled_particles.begin(), sampled_particles.end());
}

void ParticleFilter::write(std::string filename) {
    // You don't need to modify this file.
    std::ofstream dataFile;
    dataFile.open(filename, std::ios::app);
    for (int i = 0; i < particles.size(); ++i) {
        dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    dataFile.close();
}