\name{AgeBand}
\alias{AgeBand}
\title{Define an age band around a target age}
\usage{
AgeBand(TargetAge, BandWidth)
}
\arguments{
\item{TargetAge}{Target age.}
\item{BandWidth}{Number of neighboring ages to include on either side.}
}
\value{
A list containing the target age, bandwidth, lower and upper ages, and the
selected age vector.
}
\description{
Constructs the centered age band used by age-specific band-selection methods.
}
